# =============================================================================
# GP2 Progression Meta-Analysis — Model Fitting (LR + CPH), VERSION 2
# Last modified: October 10th, 2026
#
# Base formula (age + sex + disease duration)
# is fit per study and meta-analyzed (per-study fit + metafor::rma REML),
# with no other terms added yet, to see the base associations before any
# stepwise candidate terms are tested. 
#
# Term selection (Steps 0-5, per-study fit + univariate metafor::rma).
# Computes TWO independent estimates per term:
#   Meta:   fit the final model per study, pool each term's per-study
#           coefficient with a univariate metafor::rma() (REML).
#   Pooled: fit the final model ONCE on all studies' individual-level data
#           pooled together, with "study" added as a plain covariate.
# Per term, Pooled is reported only if BOTH hold (otherwise Meta is used):
#   (1) the Pooled point estimate falls inside the Meta 95% CI
#   (2) Pooled and Meta reach the same conclusion (both significant with the
#       same direction, OR both non-significant)
# Meta's tau2/I2/Cochran's Q are reported for every term regardless of which
# estimate is reported, for transparency.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(metafor)
  library(yaml)
})

# ── Load config ───────────────────────────────────────────────────────────────
config_path <- if (file.exists("config.yaml")) "config.yaml" else path.expand("~/config.yaml")
cfg <- read_yaml(config_path)

config_dir <- dirname(config_path)
resolve_cfg_path <- function(p) if (grepl("^(/|~)", p)) path.expand(p) else file.path(config_dir, p)
cfg$lr_input_file        <- resolve_cfg_path(cfg$lr_input_file)
cfg$cph_input_file       <- resolve_cfg_path(cfg$cph_input_file)
cfg$progression_log_dir  <- resolve_cfg_path(cfg$progression_log_dir)
cfg$progression_temp_dir <- resolve_cfg_path(cfg$progression_temp_dir)

# =============================================================================
# LOGGING
# =============================================================================

LOG_DIR  <- cfg$progression_log_dir
TEMP_DIR <- cfg$progression_temp_dir
dir.create(LOG_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(TEMP_DIR, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(LOG_DIR, sprintf("05_fit_progression_models_%s.log", format(Sys.time(), "%Y-%m-%d_%H%M%S")))
log_con  <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
cat(sprintf("Log started %s -> %s\n\n", format(Sys.time()), basename(log_file)))
on.exit({ sink(); close(log_con) }, add = TRUE)

# =============================================================================
# GLOBALS
# =============================================================================

ALPHA                  <- 0.05
SIGNIFICANCE_THRESHOLD <- 0.05
TODAY                  <- format(Sys.Date(), "%Y-%m-%d")

TERM_AGE_SQ   <- "age_c_sq"
TERM_DUR_SQ   <- "disease_duration_c_sq"
TERM_AGE_CUBE <- "age_c_cube"
TERM_DUR_CUBE <- "disease_duration_c_cube"

# =============================================================================
# MODEL HELPERS (needed for per-study stepwise screening)
# =============================================================================

fit_glm <- function(formula_str, data) {
  tryCatch(
    glm(as.formula(formula_str), data = data, family = binomial(link = "logit")),
    error   = function(e) { message("    glm failed (", formula_str, "): ", e$message); NULL },
    warning = function(w) suppressWarnings(
      glm(as.formula(formula_str), data = data, family = binomial(link = "logit"))
    )
  )
}

log1pexp <- function(x) ifelse(x > 30, x, log1p(exp(x)))

fit_glm_firth <- function(formula_str, data) {
  tryCatch({
    fo <- as.formula(formula_str)
    mf <- model.frame(fo, data)
    y  <- as.numeric(model.response(mf))
    X  <- model.matrix(fo, mf)

    start <- tryCatch(
      coef(suppressWarnings(glm.fit(X, y, family = binomial()))),
      error = function(e) rep(0, ncol(X))
    )
    start[!is.finite(start)] <- 0

    penalized_negloglik <- function(beta) {
      eta <- as.numeric(X %*% beta)
      p   <- plogis(eta)
      w   <- pmax(p * (1 - p), 1e-10)
      info <- crossprod(X * sqrt(w))
      ll      <- sum(y * eta) - sum(log1pexp(eta))
      penalty <- 0.5 * as.numeric(determinant(info, logarithm = TRUE)$modulus)
      -(ll + penalty)
    }

    fit_opt <- optim(start, penalized_negloglik, method = "BFGS", hessian = TRUE,
                      control = list(maxit = 1000, reltol = 1e-10))
    if (fit_opt$convergence != 0) {
      message("    Firth glm did not converge cleanly (", formula_str, ")")
    }

    beta_hat <- setNames(fit_opt$par, colnames(X))
    vcov_hat <- tryCatch(solve(fit_opt$hessian), error = function(e) NULL)
    if (is.null(vcov_hat)) {
      message("    Firth glm: Hessian not invertible (", formula_str, ")")
      return(NULL)
    }

    vcov_hat <- (vcov_hat + t(vcov_hat)) / 2
    dimnames(vcov_hat) <- list(colnames(X), colnames(X))

    structure(list(coefficients = beta_hat, vcov = vcov_hat, formula = fo),
              class = "firth_fit")
  }, error = function(e) { message("    Firth glm failed (", formula_str, "): ", e$message); NULL })
}

coef.firth_fit <- function(object, ...) object$coefficients
vcov.firth_fit <- function(object, ...) object$vcov

fit_cox <- function(formula_str, data) {
  tryCatch(
    survival::coxph(as.formula(formula_str), data = data, ties = "efron"),
    error   = function(e) { message("    coxph failed (", formula_str, "): ", e$message); NULL },
    warning = function(w) suppressWarnings(
      survival::coxph(as.formula(formula_str), data = data, ties = "efron")
    )
  )
}

# Resolving interaction term naming convention
resolve_term <- function(fit, term) {
  cn <- names(coef(fit))
  if (term %in% cn) return(term)
  if (grepl(":", term, fixed = TRUE)) {
    parts   <- strsplit(term, ":", fixed = TRUE)[[1]]
    swapped <- paste(rev(parts), collapse = ":")
    if (swapped %in% cn) return(swapped)
  }
  NA_character_
}

extract_coef_se <- function(fit, term) {
  if (is.null(fit)) return(NULL)
  rt <- resolve_term(fit, term)
  if (is.na(rt)) return(NULL)
  list(beta = as.numeric(coef(fit)[rt]), se = as.numeric(sqrt(vcov(fit)[rt, rt])))
}

# Per-study effect + CI row: feeds both the Meta side of the final-model
# comparison (per-study betas pooled via run_rma()) and the forest-plot
# per-study rows in 06_plot_forest_comparison.R. Not used by the Pooled side,
# which fits once directly on individual-level pooled data instead.
extract_full_term <- function(fit, data, term, effect_label) {
  if (is.null(fit)) return(NULL)
  rt <- resolve_term(fit, term)
  if (is.na(rt)) return(NULL)

  cf <- coef(fit)
  vc <- vcov(fit)
  beta <- as.numeric(cf[rt])
  se   <- as.numeric(sqrt(vc[rt, rt]))
  if (!is.finite(beta) || !is.finite(se) || se <= 0) return(NULL)

  z_stat <- beta / se
  p_val  <- 2 * pnorm(-abs(z_stat))
  z_crit <- qnorm(0.975)

  tibble(
    term         = term,
    coef         = beta,
    se           = se,
    effect       = exp(beta),
    CI_lower     = exp(beta - z_crit * se),
    CI_upper     = exp(beta + z_crit * se),
    p_value      = p_val,
    N            = as.integer(nrow(data)),
    N_events     = as.integer(sum(data$event, na.rm = TRUE)),
    effect_label = effect_label
  )
}

run_rma <- function(betas, ses, term, alpha = ALPHA) {
  ok    <- !is.na(betas) & !is.na(ses) & ses > 0
  betas <- betas[ok]; ses <- ses[ok]; k <- length(betas)
  if (k < 2) { message(sprintf("  [%s] <2 valid studies — skipping", term)); return(NULL) }

  fit_re <- tryCatch(metafor::rma(yi = betas, sei = ses, method = "REML"),
                     error = function(e) { message("  rma() failed: ", e$message); NULL })
  if (is.null(fit_re)) return(NULL)

  z_crit <- qnorm(1 - alpha / 2)
  tibble(
    term         = term, k_studies = as.integer(k),
    beta_RE      = as.numeric(fit_re$beta),
    se_RE        = as.numeric(fit_re$se),
    p_RE         = as.numeric(fit_re$pval),
    effect_RE    = exp(as.numeric(fit_re$beta)),
    effect_RE_lo = exp(as.numeric(fit_re$beta) - z_crit * as.numeric(fit_re$se)),
    effect_RE_hi = exp(as.numeric(fit_re$beta) + z_crit * as.numeric(fit_re$se)),
    tau2         = fit_re$tau2, I2 = fit_re$I2,
    Q            = fit_re$QE,   df_Q = as.integer(k - 1),
    p_Q          = fit_re$QEp
  )
}

build_formula <- function(base_formula, extra_terms) {
  if (length(extra_terms) == 0) return(base_formula)
  paste(base_formula, "+", paste(extra_terms, collapse = " + "))
}

# =============================================================================
# MODELS  (LR: glm/binomial on cross-sectional data; CPH: coxph on
# longitudinal survival data).
# =============================================================================

lr_sex_prep <- function(sdf) {
  sdf$sex <- factor(sdf$sex, levels = c("Female", "Male"))  # already "Male"/"Female" strings
  sdf
}

cph_sex_prep <- function(sdf) {
  sdf$biological_sex_for_qc <- factor(
    ifelse(sdf$sex == 1, "Male", ifelse(sdf$sex == 2, "Female", NA)),
    levels = c("Female", "Male")
  )
  sdf
}

lr_model <- list(
  name         = "lr",
  input_file   = cfg$lr_input_file,
  study_col    = "cohort",
  base_cols    = c("event", "age_at_baseline", "sex", "disease_duration"),
  base_formula = "event ~ age_c + sex + disease_duration_c",
  sex_var      = "sex",
  sex_coef     = "sexMale",
  sex_prep     = lr_sex_prep,
  fit_fn       = fit_glm_firth,
  # Plain glm (no Firth penalty) for the Pooled fit. Firth to
  # correct small-study/separation bias per cohort; pooling ~10k+ rows
  # across cohorts makes that far less of a concern, and the feedback this
  # was built from specifies plain glm() for the pooled model.
  pooled_fit_fn = fit_glm,
  min_n        = 20,
  min_events   = 20,
  effect_label = "OR",
  # For the LR (cross-sectional) stage, event is defined at the same visit
  # hy1_baseline describes, so hy1_baseline == 1 implies event == 0.
  # Not a meaningful covariate test there, so this model is
  # CPH-only (see run_hy1_baseline_model()).
  run_hy1_baseline = FALSE
)

cph_model <- list(
  name         = "cph",
  input_file   = cfg$cph_input_file,
  study_col    = "cohort",
  base_cols    = c("event", "study_days", "age_at_baseline", "sex", "disease_duration"),
  base_formula = "Surv(study_days, event) ~ age_c + biological_sex_for_qc + disease_duration_c",
  sex_var      = "biological_sex_for_qc",
  sex_coef     = "biological_sex_for_qcMale",
  sex_prep     = cph_sex_prep,
  fit_fn       = fit_cox,
  pooled_fit_fn = fit_cox,  # no Firth variant used for CPH in this pipeline anyway
  min_n        = 20,
  min_events   = 20,
  effect_label = "HR",
  run_hy1_baseline = TRUE
)

# Term mapping to help with plotting script
canonical_key_map <- function(stage) {
  m <- c(
    "age_c"                 = "age",
    "disease_duration_c"    = "disease_duration",
    "age_c_sq"               = "age_sq",
    "disease_duration_c_sq"  = "disease_duration_sq",
    "age_c_cube"             = "age_cube",
    "disease_duration_c_cube"= "disease_duration_cube",
    "hy1_baseline"           = "hy1_baseline"
  )
  m[[stage$sex_coef]]                                       <- "sex"
  m[[paste0("age_c:", stage$sex_coef)]]                     <- "age_x_sex"
  m[[paste0("disease_duration_c:", stage$sex_coef)]]        <- "disease_duration_x_sex"
  m
}

# =============================================================================
# PER-STUDY DATA PREP 
# =============================================================================

prepare_study_data <- function(stage, study_id, sdf, grand_mean_age, grand_mean_dd) {
  sdf <- sdf[complete.cases(sdf[, stage$base_cols]), ]

  if (nrow(sdf) < stage$min_n) {
    cat(sprintf("  Skipping %s: n=%d\n", study_id, nrow(sdf)))
    return(NULL)
  }
  n_ev <- sum(sdf$event)
  if (n_ev < stage$min_events || (nrow(sdf) - n_ev) < stage$min_events) {
    cat(sprintf("  Skipping %s: events=%d non-events=%d\n",
                study_id, as.integer(n_ev), as.integer(nrow(sdf) - n_ev)))
    return(NULL)
  }

  sdf <- stage$sex_prep(sdf)
  sdf$age_c                   <- sdf$age_at_baseline - grand_mean_age
  sdf$disease_duration_c      <- sdf$disease_duration - grand_mean_dd
  sdf$age_c_sq                <- sdf$age_c ^ 2
  sdf$disease_duration_c_sq   <- sdf$disease_duration_c ^ 2
  sdf$age_c_cube               <- sdf$age_c ^ 3
  sdf$disease_duration_c_cube  <- sdf$disease_duration_c ^ 3

  sdf
}

# =============================================================================
# BASE MODEL: per-study fit + meta-analysis, before any stepwise terms
# =============================================================================
# Fits stage$base_formula (age_c + sex + disease_duration_c) per study, then
# pools each term's per-study coefficient with a univariate metafor::rma()
# (REML). Same two-stage idea as the stepwise screening steps below, but
# on the base terms jointly rather than one candidate term at a time, and
# with no candidate terms added yet. This is what the main forest-plot
# figure (06_plot_forest_comparison.R) shows for age/sex/disease duration.

run_base_model <- function(stage, study_data, valid_studies, key_map, alpha = ALPHA) {
  k <- length(valid_studies)
  base_terms <- c("age_c", stage$sex_coef, "disease_duration_c")

  cat(sprintf("\n  Base model -- per-study fit + meta-analysis: %d studies\n", k))
  cat(sprintf("  Model: %s\n", stage$base_formula))

  base_fits <- setNames(lapply(valid_studies, function(sid) stage$fit_fn(stage$base_formula, study_data[[sid]])), valid_studies)

  base_rows <- list()
  for (sid in valid_studies) {
    for (term in base_terms) {
      row <- extract_full_term(base_fits[[sid]], study_data[[sid]], term, stage$effect_label)
      if (!is.null(row)) base_rows <- c(base_rows, list(cbind(study = sid, row)))
    }
  }
  per_study_df <- bind_rows(base_rows)
  per_study_df$term_key <- key_map[per_study_df$term]

  meta_rows <- lapply(base_terms, function(term) {
    ests <- lapply(valid_studies, function(sid) extract_coef_se(base_fits[[sid]], term))
    run_rma(sapply(ests, function(x) if (is.null(x)) NA_real_ else x$beta),
            sapply(ests, function(x) if (is.null(x)) NA_real_ else x$se),
            term, alpha = alpha)
  })
  base_df <- bind_rows(meta_rows)
  base_df$term_key <- key_map[base_df$term]

  list(base_df = base_df, per_study_df = per_study_df)
}

# =============================================================================
# HY1_BASELINE MODEL: base + hy1_baseline (CPH only; see run_hy1_baseline
# flag on cph_model/lr_model above)
# =============================================================================
# A separate, fixed model (not part of the Step 0-5 stepwise selection),
# same two-stage per-study-fit + univariate rma() methodology as
# run_base_model() above -- just with hy1_baseline added to the formula.

run_hy1_baseline_model <- function(stage, study_data, valid_studies, key_map, alpha = ALPHA) {
  hy1_formula <- paste(stage$base_formula, "+ hy1_baseline")
  hy1_terms   <- c("age_c", stage$sex_coef, "disease_duration_c", "hy1_baseline")

  cat(sprintf("\n  HY1-baseline model -- per-study fit + meta-analysis\n"))
  cat(sprintf("  Model: %s\n", hy1_formula))

  # Restrict to studies with complete hy1_baseline data, re-checking the
  # inclusion floor (n/events/non-events) after that restriction. A study
  # that passed prepare_study_data()'s check on the base columns can still
  # drop below the floor once rows missing hy1_baseline are excluded.
  hy1_study_data <- list()
  for (sid in valid_studies) {
    sdf <- study_data[[sid]]
    sdf <- sdf[complete.cases(sdf[, "hy1_baseline", drop = FALSE]), ]
    n_ev <- sum(sdf$event)
    if (nrow(sdf) < stage$min_n || n_ev < stage$min_events || (nrow(sdf) - n_ev) < stage$min_events) {
      cat(sprintf("  Skipping %s: n=%d events=%d non-events=%d after restricting to complete hy1_baseline\n",
                  sid, nrow(sdf), as.integer(n_ev), as.integer(nrow(sdf) - n_ev)))
      next
    }
    hy1_study_data[[sid]] <- sdf
  }
  hy1_valid_studies <- names(hy1_study_data)
  cat(sprintf("  %d/%d studies retained (hy1_baseline available and inclusion floor still met).\n",
              length(hy1_valid_studies), length(valid_studies)))

  if (length(hy1_valid_studies) < 2) {
    cat(sprintf("[%s] HY1-baseline model: fewer than 2 studies available -- skipping.\n", stage$name))
    return(NULL)
  }

  hy1_fits <- setNames(lapply(hy1_valid_studies, function(sid) stage$fit_fn(hy1_formula, hy1_study_data[[sid]])), hy1_valid_studies)

  hy1_rows <- list()
  for (sid in hy1_valid_studies) {
    for (term in hy1_terms) {
      row <- extract_full_term(hy1_fits[[sid]], hy1_study_data[[sid]], term, stage$effect_label)
      if (!is.null(row)) hy1_rows <- c(hy1_rows, list(cbind(study = sid, row)))
    }
  }
  per_study_df <- bind_rows(hy1_rows)
  per_study_df$term_key <- key_map[per_study_df$term]

  meta_rows <- lapply(hy1_terms, function(term) {
    ests <- lapply(hy1_valid_studies, function(sid) extract_coef_se(hy1_fits[[sid]], term))
    run_rma(sapply(ests, function(x) if (is.null(x)) NA_real_ else x$beta),
            sapply(ests, function(x) if (is.null(x)) NA_real_ else x$se),
            term, alpha = alpha)
  })
  hy1_df <- bind_rows(meta_rows)
  hy1_df$term_key <- key_map[hy1_df$term]

  list(hy1_df = hy1_df, per_study_df = per_study_df)
}

# =============================================================================
# FINAL MODEL: per-term Pooled-vs-Meta comparison
# =============================================================================
# For every term in the final model (base terms + whatever stepwise
# screening confirmed), two independent estimates are computed and compared

run_final_comparison <- function(stage, study_data, valid_studies, confirmed_terms, confirmed_coefs, key_map, alpha = ALPHA) {
  k <- length(valid_studies)
  final_formula     <- build_formula(stage$base_formula, confirmed_terms)
  joint_terms_final <- c("age_c", stage$sex_coef, "disease_duration_c", confirmed_coefs)

  cat(sprintf("\n  Final model -- Pooled vs. Meta comparison: %d studies\n", k))
  cat(sprintf("  Model: %s\n", final_formula))

  # --- Per-study fits of the final model: feed the Meta side (per-study
  # betas pooled with run_rma()) and the forest-plot per-study rows. ---
  final_fits <- setNames(lapply(valid_studies, function(sid) stage$fit_fn(final_formula, study_data[[sid]])), valid_studies)

  final_rows <- list()
  for (sid in valid_studies) {
    for (term in joint_terms_final) {
      row <- extract_full_term(final_fits[[sid]], study_data[[sid]], term, stage$effect_label)
      if (!is.null(row)) final_rows <- c(final_rows, list(cbind(study = sid, row)))
    }
  }
  per_study_df <- bind_rows(final_rows)

  meta_list <- setNames(lapply(joint_terms_final, function(term) {
    ests <- lapply(valid_studies, function(sid) extract_coef_se(final_fits[[sid]], term))
    run_rma(sapply(ests, function(x) if (is.null(x)) NA_real_ else x$beta),
            sapply(ests, function(x) if (is.null(x)) NA_real_ else x$se),
            term, alpha = alpha)
  }), joint_terms_final)

  # --- Pooled: one fit directly on all studies' individual-level data,
  # pooled together, with "study" added as a covariate. ---
  pooled_df <- bind_rows(lapply(valid_studies, function(sid) {
    sdf <- study_data[[sid]]
    sdf$study <- sid
    sdf
  }))
  pooled_df$study <- factor(pooled_df$study)

  pooled_formula_str <- paste(final_formula, "+ study")
  cat(sprintf("  Pooled model: %s\n", pooled_formula_str))

  model_vars     <- all.vars(as.formula(pooled_formula_str))
  missing_counts <- sapply(pooled_df[, model_vars], function(x) sum(is.na(x)))
  if (any(missing_counts > 0)) {
    cat("  Dropping rows with NA in model variables before pooled fit:\n")
    print(missing_counts[missing_counts > 0])
  }
  n_before  <- nrow(pooled_df)
  pooled_df <- pooled_df[complete.cases(pooled_df[, model_vars]), ]
  if (nrow(pooled_df) < n_before) {
    cat(sprintf("  %d/%d rows retained after dropping incomplete cases\n", nrow(pooled_df), n_before))
  }

  pooled_fit <- stage$pooled_fit_fn(pooled_formula_str, pooled_df)
  if (is.null(pooled_fit)) {
    message(sprintf("  [%s] Pooled model fit failed -- every term will fall back to Meta.", stage$name))
  }

  z_crit <- qnorm(1 - alpha / 2)

  rows <- lapply(joint_terms_final, function(term) {
    meta <- meta_list[[term]]
    if (is.null(meta)) {
      return(tibble(
        term = term, k_studies = NA_integer_,
        meta_beta = NA_real_, meta_se = NA_real_, meta_p = NA_real_,
        meta_effect = NA_real_, meta_effect_lo = NA_real_, meta_effect_hi = NA_real_,
        meta_tau2 = NA_real_, meta_I2 = NA_real_, meta_Qp = NA_real_,
        pooled_beta = NA_real_, pooled_se = NA_real_, pooled_p = NA_real_,
        pooled_effect = NA_real_, pooled_effect_lo = NA_real_, pooled_effect_hi = NA_real_,
        ci_overlap = NA, same_conclusion = NA, decision = NA_character_,
        beta_RE = NA_real_, se_RE = NA_real_, p_RE = NA_real_,
        effect_RE = NA_real_, effect_RE_lo = NA_real_, effect_RE_hi = NA_real_
      ))
    }

    pest <- if (is.null(pooled_fit)) NULL else extract_coef_se(pooled_fit, term)

    if (is.null(pest)) {
      # Pooled unavailable (fit failed, or term not resolvable in it) --
      # nothing to compare against, so Meta is reported by definition.
      decision        <- "Meta"
      pooled_beta     <- NA_real_; pooled_se <- NA_real_; pooled_p <- NA_real_
      ci_overlap      <- NA
      same_conclusion <- NA
    } else {
      pooled_beta <- pest$beta
      pooled_se   <- pest$se
      pooled_p    <- 2 * pnorm(-abs(pooled_beta / pooled_se))

      ci_overlap <- pooled_beta >= (meta$beta_RE - z_crit * meta$se_RE) &&
                    pooled_beta <= (meta$beta_RE + z_crit * meta$se_RE)

      meta_sig   <- meta$p_RE < SIGNIFICANCE_THRESHOLD
      pooled_sig <- pooled_p  < SIGNIFICANCE_THRESHOLD
      same_dir   <- sign(meta$beta_RE) == sign(pooled_beta)
      same_conclusion <- (meta_sig && pooled_sig && same_dir) || (!meta_sig && !pooled_sig)

      decision <- if (isTRUE(ci_overlap) && isTRUE(same_conclusion)) "Pooled" else "Meta"
    }

    reported_beta <- if (decision == "Pooled") pooled_beta else meta$beta_RE
    reported_se   <- if (decision == "Pooled") pooled_se   else meta$se_RE
    reported_p    <- if (decision == "Pooled") pooled_p    else meta$p_RE

    tibble(
      term             = term,
      k_studies        = meta$k_studies,
      meta_beta        = meta$beta_RE,   meta_se = meta$se_RE,   meta_p = meta$p_RE,
      meta_effect      = meta$effect_RE, meta_effect_lo = meta$effect_RE_lo, meta_effect_hi = meta$effect_RE_hi,
      meta_tau2        = meta$tau2,      meta_I2 = meta$I2,      meta_Qp = meta$p_Q,
      pooled_beta      = pooled_beta,    pooled_se = pooled_se,  pooled_p = pooled_p,
      pooled_effect    = if (is.na(pooled_beta)) NA_real_ else exp(pooled_beta),
      pooled_effect_lo = if (is.na(pooled_beta)) NA_real_ else exp(pooled_beta - z_crit * pooled_se),
      pooled_effect_hi = if (is.na(pooled_beta)) NA_real_ else exp(pooled_beta + z_crit * pooled_se),
      ci_overlap       = ci_overlap,
      same_conclusion  = same_conclusion,
      decision         = decision,
      beta_RE          = reported_beta,
      se_RE            = reported_se,
      p_RE             = reported_p,
      effect_RE        = exp(reported_beta),
      effect_RE_lo     = exp(reported_beta - z_crit * reported_se),
      effect_RE_hi     = exp(reported_beta + z_crit * reported_se)
    )
  })

  final_df <- bind_rows(rows)
  final_df$term_key <- key_map[final_df$term]

  list(final_df = final_df, per_study_df = per_study_df)
}

# =============================================================================
# MAIN PIPELINE  
# =============================================================================

run_pipeline <- function(stage) {
  cat(sprintf("\n\n%s\n%s PIPELINE (one-stage comparison)\n%s\n", strrep("#", 70), toupper(stage$name), strrep("#", 70)))

  df <- read.csv(stage$input_file, stringsAsFactors = FALSE)
  df <- df[df[[stage$study_col]] != "TRD", ]
  cat(sprintf("Loaded %s\n  %d rows, %d cohorts (TRD excluded)\n",
              basename(stage$input_file), nrow(df), length(unique(df[[stage$study_col]]))))

  grand_mean_age <- mean(df$age_at_baseline,  na.rm = TRUE)
  grand_mean_dd  <- mean(df$disease_duration, na.rm = TRUE)
  cat(sprintf("Grand mean age at baseline: %.4f | Grand mean disease duration: %.4f\n\n",
              grand_mean_age, grand_mean_dd))

  studies    <- unique(df[[stage$study_col]])
  study_data <- list()
  for (sid in studies) {
    cat(sprintf("%s\nStudy: %s\n", strrep("=", 60), sid))
    sdf <- prepare_study_data(stage, sid, df[df[[stage$study_col]] == sid, ], grand_mean_age, grand_mean_dd)
    if (!is.null(sdf)) study_data[[sid]] <- sdf
  }
  valid_studies <- names(study_data)
  cat(sprintf("\n%s\n%d/%d studies passed inclusion (n>=%d, events & non-events>=%d).\n\n",
              strrep("=", 60), length(valid_studies), length(studies), stage$min_n, stage$min_events))

  key_map <- canonical_key_map(stage)

  # --- Base model: per-study fit + meta-analysis of age/sex/disease duration
  # alone, before any stepwise candidate terms are tested. Feeds the main
  # forest-plot figure. ---
  base_result <- run_base_model(stage, study_data, valid_studies, key_map)
  base_df     <- base_result$base_df

  cat("\nBase model -- meta-analysis (pre-stepwise):\n")
  print(base_df[, c("term", "k_studies", "beta_RE", "se_RE", "p_RE",
                     "effect_RE", "effect_RE_lo", "effect_RE_hi", "I2", "p_Q")])

  write.csv(base_df,
            file.path(TEMP_DIR, sprintf("%s_base_model_%s.csv", stage$name, TODAY)), row.names = FALSE)
  cat(sprintf("[%s] Base-model output saved to %s/\n", stage$name, basename(TEMP_DIR)))

  write.csv(base_result$per_study_df,
            file.path(TEMP_DIR, sprintf("%s_base_model_per_study_%s.csv", stage$name, TODAY)), row.names = FALSE)
  cat(sprintf("[%s] Base-model per-study output saved to %s/\n", stage$name, basename(TEMP_DIR)))

  # --- HY1-baseline model: base + hy1_baseline (CPH only) ---------------------
  if (isTRUE(stage$run_hy1_baseline)) {
    hy1_result <- run_hy1_baseline_model(stage, study_data, valid_studies, key_map)
    if (!is.null(hy1_result)) {
      cat("\nHY1-baseline model -- meta-analysis:\n")
      print(hy1_result$hy1_df[, c("term", "k_studies", "beta_RE", "se_RE", "p_RE",
                                   "effect_RE", "effect_RE_lo", "effect_RE_hi", "I2", "p_Q")])

      write.csv(hy1_result$hy1_df,
                file.path(TEMP_DIR, sprintf("%s_hy1_baseline_model_%s.csv", stage$name, TODAY)), row.names = FALSE)
      write.csv(hy1_result$per_study_df,
                file.path(TEMP_DIR, sprintf("%s_hy1_baseline_model_per_study_%s.csv", stage$name, TODAY)), row.names = FALSE)
      cat(sprintf("[%s] HY1-baseline model outputs saved to %s/\n", stage$name, basename(TEMP_DIR)))
    }
  }

  # --- Stepwise term selection
  confirmed_terms <- character(0)
  confirmed_coefs <- character(0)
  sig_flags       <- list()

  # Every step's univariate pooled result captured
  # regardless of significance and written to a separate decision-steps CSV so
  # the sex interactions (and any other non-selected term) are still visible
  decision_rows <- list()

  # Per-study effect + CI for every candidate term tested, by
  # coef_name; used only to give the
  # forest plot per-study rows for terms that don't end up in the final
  # model (see extract_full_term()).
  all_candidate_per_study <- list()

  run_step <- function(step_key, step_label, formula_frag, coef_name, eligible = TRUE) {
    step_formula <- build_formula(stage$base_formula, confirmed_terms)
    test_formula <- paste(step_formula, "+", formula_frag)
    cat(sprintf("\n%s\n%s (%s)\n%s\n", strrep("-", 60), step_label, coef_name, strrep("-", 60)))
    cat(sprintf("  Model: %s\n", test_formula))

    fits <- setNames(lapply(valid_studies, function(sid) stage$fit_fn(test_formula, study_data[[sid]])), valid_studies)
    ests <- lapply(valid_studies, function(sid) extract_coef_se(fits[[sid]], coef_name))
    meta <- run_rma(sapply(ests, function(x) if (is.null(x)) NA_real_ else x$beta),
                     sapply(ests, function(x) if (is.null(x)) NA_real_ else x$se),
                     coef_name)
    if (!is.null(meta)) print(meta[, c("term", "k_studies", "beta_RE", "se_RE", "p_RE", "I2", "p_Q")])

    cand_rows <- list()
    for (sid in valid_studies) {
      row <- extract_full_term(fits[[sid]], study_data[[sid]], coef_name, stage$effect_label)
      if (!is.null(row)) cand_rows <- c(cand_rows, list(cbind(study = sid, row)))
    }
    if (length(cand_rows) > 0) all_candidate_per_study[[coef_name]] <<- bind_rows(cand_rows)

    sig          <- !is.null(meta) && meta$p_RE < SIGNIFICANCE_THRESHOLD
    add_to_model <- sig && eligible
    if (!eligible) {
      cat(sprintf("  %s: p=%s -- quadratic counterpart not significant, do not include; reporting estimate only\n",
                  coef_name, if (!is.null(meta)) sprintf("%.4f", meta$p_RE) else "NA"))
    } else {
      cat(sprintf("  %s %s (p=%.4f)\n", coef_name,
                  if (sig) "significant -> added to model" else "NS -> not added",
                  if (!is.null(meta)) meta$p_RE else NA))
    }

    if (add_to_model) {
      confirmed_terms <<- c(confirmed_terms, formula_frag)
      confirmed_coefs <<- c(confirmed_coefs, coef_name)
    }
    sig_flags[[step_key]] <<- sig
    decision_rows[[step_label]] <<- if (!is.null(meta)) mutate(meta, step = step_label) else NULL
    invisible(add_to_model)
  }

  run_step("dur_sq",  "Step0_disease_duration_c_sq", TERM_DUR_SQ, TERM_DUR_SQ)
  run_step("age_sq",  "Step1_age_c_sq",              TERM_AGE_SQ, TERM_AGE_SQ)
  run_step("age_sex", "Step2_age_c_x_sex",
           paste0("age_c:", stage$sex_var), paste0("age_c:", stage$sex_coef))
  run_step("dur_sex", "Step3_disease_duration_c_x_sex",
           paste0("disease_duration_c:", stage$sex_var), paste0("disease_duration_c:", stage$sex_coef))
  run_step("dur_cube", "Step4_disease_duration_c_cube", TERM_DUR_CUBE, TERM_DUR_CUBE,
           eligible = isTRUE(sig_flags[["dur_sq"]]))
  run_step("age_cube", "Step5_age_c_cube", TERM_AGE_CUBE, TERM_AGE_CUBE,
           eligible = isTRUE(sig_flags[["age_sq"]]))

  cat(sprintf("\n%s\nFinal model: %s\n%s\n\n", strrep("*", 60),
              build_formula(stage$base_formula, confirmed_terms), strrep("*", 60)))

  # --- Final step: per-term Pooled-vs-Meta comparison ---
  final_formula     <- build_formula(stage$base_formula, confirmed_terms)
  joint_terms_final <- c("age_c", stage$sex_coef, "disease_duration_c", confirmed_coefs)

  final_result <- run_final_comparison(stage, study_data, valid_studies, confirmed_terms, confirmed_coefs, key_map)
  final_df     <- final_result$final_df
  final_df$in_final_model <- TRUE

  cat("\nFinal model -- Pooled vs. Meta comparison:\n")
  print(final_df[, c("term", "k_studies", "meta_p", "pooled_p", "ci_overlap", "same_conclusion", "decision",
                      "effect_RE", "effect_RE_lo", "effect_RE_hi", "p_RE")])

  write.csv(final_df,
            file.path(TEMP_DIR, sprintf("%s_onestage_joint_%s.csv", stage$name, TODAY)), row.names = FALSE)
  cat(sprintf("\n[%s] Final comparison output saved to %s/\n", stage$name, basename(TEMP_DIR)))

  # --- Per-study rows for forest-plot display (from run_final_comparison()'s
  # per-study final-model fits -- 
  per_study_df <- final_result$per_study_df
  per_study_df$term_key       <- key_map[per_study_df$term]
  per_study_df$in_final_model <- TRUE

  # Also include per-study rows for candidate terms NOT selected into the
  # final model 
  non_final_coefs <- setdiff(names(all_candidate_per_study), confirmed_coefs)
  if (length(non_final_coefs) > 0) {
    non_final_study_df <- bind_rows(lapply(non_final_coefs, function(cn) {
      sdf <- all_candidate_per_study[[cn]]
      if (is.null(sdf)) return(NULL)
      sdf$term_key       <- key_map[sdf$term]
      sdf$in_final_model <- FALSE
      sdf
    }))
    if (nrow(non_final_study_df) > 0) per_study_df <- bind_rows(per_study_df, non_final_study_df)
  }

  write.csv(per_study_df,
            file.path(TEMP_DIR, sprintf("%s_onestage_per_study_%s.csv", stage$name, TODAY)), row.names = FALSE)
  cat(sprintf("[%s] Per-study output saved to %s/\n", stage$name, basename(TEMP_DIR)))

  # --- Decision-steps log: every candidate term's univariate pooled result,
  # significant or not, kept separate from the
  # final model comparison output above.
  decision_summary <- bind_rows(decision_rows)
  if (nrow(decision_summary) > 0) decision_summary$term_key <- key_map[decision_summary$term]
  write.csv(decision_summary,
            file.path(TEMP_DIR, sprintf("%s_decision_steps_%s.csv", stage$name, TODAY)), row.names = FALSE)
  cat(sprintf("[%s] Decision-steps log saved to %s/\n", stage$name, basename(TEMP_DIR)))

  invisible(list(
    stage_name         = stage$name,
    final_formula      = final_formula,
    joint_terms_final  = joint_terms_final,
    sig_flags          = sig_flags,
    base_model         = base_df,
    final_comparison   = final_df,
    per_study          = per_study_df,
    decision_steps     = decision_summary,
    grand_mean_age     = grand_mean_age,
    grand_mean_dd      = grand_mean_dd
  ))
}

# =============================================================================
# COMBINED SUMMARY TABLE
# =============================================================================
# The analysis is (0) base-model meta-analysis (pre-stepwise), (1) stepwise
# term selection, (2) meta-analyze the final model (per-study fit + rma()),
# (3) pool data and re-fit the final formula once. This table is where steps
# (2) and (3) are compared side by side, one row per final-model term per
# stage. 

build_summary_table <- function(stage_name, effect_label, final_df) {
  tibble(
    model         = stage_name,
    term          = final_df$term_key,
    effect_label  = effect_label,   # "OR" for lr rows, "HR" for cph rows
    metaRE_beta   = final_df$meta_beta,
    metaRE_se     = final_df$meta_se,
    metaRE_effect = final_df$meta_effect,
    metaRE_CI     = sprintf("%.3f–%.3f", final_df$meta_effect_lo, final_df$meta_effect_hi),
    I2            = final_df$meta_I2,
    Qp            = final_df$meta_Qp,
    Pooled_beta   = final_df$pooled_beta,
    Pooled_SE     = final_df$pooled_se,
    Pooled_effect = final_df$pooled_effect,
    Pooled_CI     = ifelse(is.na(final_df$pooled_beta), NA_character_,
                            sprintf("%.3f–%.3f", final_df$pooled_effect_lo, final_df$pooled_effect_hi))
  )
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  lr_out  <- run_pipeline(lr_model)
  cph_out <- run_pipeline(cph_model)

  cat(sprintf("\n\n%s\nSUMMARY (v2 -- Pooled-vs-Meta comparison)\n%s\n", strrep("#", 70), strrep("#", 70)))
  cat(sprintf("[LR]  final model: %s\n", lr_out$final_formula))
  cat(sprintf("[CPH] final model: %s\n", cph_out$final_formula))
  cat("\nPer-term decision (Pooled vs. Meta) is in the \"decision\" column of\n")
  cat("<stage>_onestage_joint_<date>.csv. Compare against <stage>_joint_meta_<date>.csv\n")
  cat("(v1, two-stage rma.mv) for the original methodology.\n")

  summary_table <- bind_rows(
    build_summary_table("lr",  lr_model$effect_label,  lr_out$final_comparison),
    build_summary_table("cph", cph_model$effect_label, cph_out$final_comparison)
  )
  summary_path <- file.path(TEMP_DIR, sprintf("progression_final_summary_%s.csv", TODAY))
  write.csv(summary_table, summary_path, row.names = FALSE)
  cat(sprintf("\nCombined Meta-vs-Pooled summary table saved to %s\n", basename(summary_path)))
  print(summary_table)
}

main()
