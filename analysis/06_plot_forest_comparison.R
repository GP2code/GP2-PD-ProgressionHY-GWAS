# =============================================================================
# GP2 Progression Meta-Analysis — LR vs. CPH Forest Plot Comparison
#
# Visualization after running 05_fit_progression_models.R.
# Per-study rows are drawn and each study's effect + CI comes from fitting the
# relevant model per study (base model for the main figure, final model for
# everything else).
#
# The diamond at the bottom of each panel is a Meta RE estimate, two-stage
# per-study fit + rma(), but which fit differs by figure:
#
#   - Main figure (age, sex, disease duration): the BASE model (age + sex +
#     disease duration only, no stepwise terms added), from
#     <stage>_base_model_<date>.csv / <stage>_base_model_per_study_<date>.csv.
#     This shows the pre-stepwise associations, before adjustment for
#     whatever quadratic/cubic/interaction terms stepwise selection added.
#   - Supplementary figures (one per stepwise-tested term): the FINAL model
#     if selected (meta_effect/meta_effect_lo/meta_effect_hi etc. from
#     <stage>_onestage_joint_<date>.csv), or the univariate screening
#     estimate if not selected (from <stage>_decision_steps_<date>.csv).
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(grDevices)
  library(yaml)
})

# ── Load config ───────────────────────────────────────────────────────────────
config_path <- if (file.exists("config.yaml")) "config.yaml" else path.expand("~/config.yaml")
cfg <- read_yaml(config_path)

config_dir <- dirname(config_path)
resolve_cfg_path <- function(p) if (grepl("^(/|~)", p)) path.expand(p) else file.path(config_dir, p)
cfg$lr_input_file           <- resolve_cfg_path(cfg$lr_input_file)
cfg$progression_log_dir     <- resolve_cfg_path(cfg$progression_log_dir)
cfg$progression_temp_dir    <- resolve_cfg_path(cfg$progression_temp_dir)
cfg$progression_figures_dir <- resolve_cfg_path(cfg$progression_figures_dir)

# =============================================================================
# LOGGING
# =============================================================================

LOG_DIR     <- cfg$progression_log_dir
TEMP_DIR    <- cfg$progression_temp_dir
FIGURES_DIR <- cfg$progression_figures_dir
dir.create(LOG_DIR,     showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(LOG_DIR, sprintf("06_plot_forest_comparison_%s.log", format(Sys.time(), "%Y-%m-%d_%H%M%S")))
log_con  <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
cat(sprintf("Log started %s -> %s\n\n", format(Sys.time()), basename(log_file)))
on.exit({ sink(); close(log_con) }, add = TRUE)

# =============================================================================
# GLOBALS
# =============================================================================

TODAY <- format(Sys.Date(), "%Y-%m-%d")

COL_RE           <- "#9C27B0"
COL_NOT_SELECTED <- "#9E9E9E"
COL_BASE         <- "#00897B"
COL_HY1          <- "#C62828"

BASE_TERM_KEYS <- c("age", "sex", "disease_duration")

# Per-study rows are colored by cohort disease-stage category
# The shorter/longer disease-duration are based on median-split cutoff from the
# data, so the category names always reflect the real split point used.
CATEGORY_COLOR_VALUES <- c("#2E8B57", "#1E88E5", "#E64A19", "#8D6E63")
CATEGORY_COLORS <- setNames(CATEGORY_COLOR_VALUES,
                            c("De novo", "Shorter disease duration", "Longer disease duration", "Other"))
CATEGORY_ORDER <- c("De novo", "Shorter disease duration", "Longer disease duration", "Other")

COHORT_ANNOTATIONS <- c(
  "IPDGCAF-NG" = "†",
  "Nigeria-PD" = "†",
  "NTUH"       = "‡"
)
COHORT_FOOTNOTE <- "† African/African Admixed cohort     ‡ East Asian cohort"

TERM_LABEL_MAP <- c(
  age                     = "Age (baseline, centered)",
  sex                     = "Sex (Male)",
  disease_duration        = "Disease Duration (centered)",
  age_sq                  = "Age²",
  disease_duration_sq     = "Disease Duration²",
  age_cube                = "Age³",
  disease_duration_cube   = "Disease Duration³",
  age_x_sex               = "Age × Sex (Male)",
  disease_duration_x_sex  = "Disease Duration × Sex (Male)",
  hy1_baseline            = "HY = 1 at Baseline"
)
format_term <- function(term_key) {
  if (term_key %in% names(TERM_LABEL_MAP)) TERM_LABEL_MAP[[term_key]] else term_key
}
format_effect_ci <- function(effect, lo, hi) {
  span     <- max(hi - lo, 1e-6)
  decimals <- max(3, min(4, ceiling(-log10(span)) + 1))
  fmt      <- sprintf("%%.%df (%%.%df–%%.%df)", decimals, decimals, decimals)
  sprintf(fmt, effect, lo, hi)
}

MODEL_SPECS <- list(
  lr  = list(label = "LR",  effect_label = "Odds Ratio",    effect_axis = "Odds Ratio (95% CI)"),
  cph = list(label = "CPH", effect_label = "Hazard Ratio",  effect_axis = "Hazard Ratio (95% CI)")
)

PANEL_WIDTH_IN       <- 22
ROW_UNIT_IN_SINGLE   <- 0.70
MARGIN_OVERHEAD_IN   <- 5.7
LEGEND_RESERVE_UNITS <- 6.5

n_term_rows <- function(studies, cohort_class_map) {
  studies <- unique(studies)
  if (length(studies) == 0) return(0)
  classes <- cohort_class_map[studies]
  classes[is.na(classes)] <- "Other"
  length(studies) + length(unique(classes))
}

# =============================================================================
# LOAD LATEST TEMP/ OUTPUTS
# =============================================================================

find_latest <- function(dir, filename_prefix) {
  files <- list.files(dir, pattern = sprintf("^%s_\\d{4}-\\d{2}-\\d{2}\\.csv$", filename_prefix))
  if (length(files) == 0) stop(sprintf("No files found matching %s_<date>.csv in %s", filename_prefix, dir))
  dates <- as.Date(sub(".*_(\\d{4}-\\d{2}-\\d{2})\\.csv$", "\\1", files))
  file.path(dir, files[which.max(dates)])
}

load_model_results <- function(name) {
  onestage_path <- find_latest(TEMP_DIR, sprintf("%s_onestage_joint", name))
  study_path    <- find_latest(TEMP_DIR, sprintf("%s_onestage_per_study", name))
  decision_path <- find_latest(TEMP_DIR, sprintf("%s_decision_steps", name))
  cat(sprintf("[%s] onestage joint:  %s\n", name, basename(onestage_path)))
  cat(sprintf("[%s] onestage study:  %s\n", name, basename(study_path)))
  cat(sprintf("[%s] decision steps:  %s\n", name, basename(decision_path)))

  onestage_df <- read.csv(onestage_path, stringsAsFactors = FALSE)
  decision_df <- read.csv(decision_path, stringsAsFactors = FALSE)

  non_final <- decision_df[!(decision_df$term_key %in% onestage_df$term_key), ]
  if (nrow(non_final) > 0) {
    non_final$pooling         <- "univariate_rma"
    non_final$in_final_model  <- FALSE
  }
  joint_df <- bind_rows(onestage_df, non_final)

  list(
    joint = joint_df,
    study = read.csv(study_path, stringsAsFactors = FALSE)
  )
}

# Base-model outputs (pre-stepwise age + sex + disease duration, per-study
# fit + rma()) 
load_base_model <- function(name) {
  base_path  <- find_latest(TEMP_DIR, sprintf("%s_base_model", name))
  study_path <- find_latest(TEMP_DIR, sprintf("%s_base_model_per_study", name))
  cat(sprintf("[%s] base model:       %s\n", name, basename(base_path)))
  cat(sprintf("[%s] base model study: %s\n", name, basename(study_path)))

  list(
    joint = read.csv(base_path, stringsAsFactors = FALSE),
    study = read.csv(study_path, stringsAsFactors = FALSE)
  )
}

# HY1-baseline model outputs (base + hy1_baseline, CPH only) 
load_hy1_baseline_model <- function(name) {
  hy1_path   <- find_latest(TEMP_DIR, sprintf("%s_hy1_baseline_model", name))
  study_path <- find_latest(TEMP_DIR, sprintf("%s_hy1_baseline_model_per_study", name))
  cat(sprintf("[%s] hy1_baseline model:       %s\n", name, basename(hy1_path)))
  cat(sprintf("[%s] hy1_baseline model study: %s\n", name, basename(study_path)))

  list(
    joint = read.csv(hy1_path, stringsAsFactors = FALSE),
    study = read.csv(study_path, stringsAsFactors = FALSE)
  )
}

# =============================================================================
# COHORT CLASSIFICATION  (de novo / shorter / longer disease duration, split
# at the median cohort-level disease duration, not a fixed year cutoff)
# =============================================================================

classify_cohorts <- function() {
  lr_cs_path <- find_latest(dirname(cfg$lr_input_file), "clin_hy_cs_dataset_for_LR")
  cat(sprintf("Cohort classification source: %s\n", basename(lr_cs_path)))

  df <- read.csv(lr_cs_path, stringsAsFactors = FALSE)
  df <- df[df$cohort != "TRD", ]

  cohort_means <- df %>%
    group_by(cohort) %>%
    summarise(mean_dd = mean(disease_duration, na.rm = TRUE), .groups = "drop")

  # Median split over non-PPMI cohorts only. PPMI is a de novo cohort
  cutoff <- median(cohort_means$mean_dd[cohort_means$cohort != "PPMI"], na.rm = TRUE)
  cat(sprintf("Disease-duration median-split cutoff (excl. PPMI): %.2f years\n", cutoff))

  lbl_lo <- sprintf("≤ %.1f yrs disease duration", cutoff)
  lbl_hi <- sprintf("> %.1f yrs disease duration", cutoff)

  CATEGORY_COLORS <<- setNames(CATEGORY_COLOR_VALUES, c("De novo", lbl_lo, lbl_hi, "Other"))
  CATEGORY_ORDER  <<- c("De novo", lbl_lo, lbl_hi, "Other")

  summary_df <- cohort_means %>%
    mutate(cohort_class = case_when(
      cohort == "PPMI"  ~ "De novo",
      mean_dd <= cutoff ~ lbl_lo,
      mean_dd >  cutoff ~ lbl_hi,
      TRUE              ~ "Other"
    ))

  cat("Cohort classification:\n")
  print(summary_df[, c("cohort", "mean_dd", "cohort_class")])

  setNames(summary_df$cohort_class, summary_df$cohort)
}

# =============================================================================
# PANEL DRAWING  (one column: per-study rows + pooled diamond)
# =============================================================================

draw_empty_panel <- function(title) {
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "", main = title, cex.main = 0.85)
  text(0.5, 0.5, "Term not selected\nin final model", cex = 1.0, col = "grey50")
}

draw_panel <- function(term_key, joint_df, study_df, spec, cohort_class_map, term_label = NULL, kind = NULL) {
  jrow  <- joint_df[joint_df$term_key == term_key, ]
  title <- if (is.null(term_label)) sprintf("%s [%s]", spec$label, spec$effect_label) else
                                     sprintf("%s\n%s [%s]", term_label, spec$label, spec$effect_label)

  if (nrow(jrow) == 0) {
    draw_empty_panel(title)
    return(invisible(NULL))
  }
  jrow <- jrow[1, ]

  # The diamond is always a Meta RE estimate (two-stage: per-study fit +
  # univariate rma()), but which fit it comes from depends on `kind`:
  #   "base"     -- pre-stepwise base model (age+sex+disease duration only),
  #                 from <stage>_base_model_<date>.csv (run_rma() shape).
  #   "final"    -- term selected into the final model, from meta_* columns
  #                 in <stage>_onestage_joint_<date>.csv.
  #   "screened" -- term tested during stepwise screening but not selected,
  #                 from <stage>_decision_steps_<date>.csv (run_rma() shape,
  #                 same columns as "base").
  #   "hy1"      -- base + hy1_baseline model (CPH only), from
  #                 <stage>_hy1_baseline_model_<date>.csv (run_rma() shape,
  #                 same columns as "base").
  if (is.null(kind)) kind <- if (isTRUE(jrow$in_final_model)) "final" else "screened"

  if (kind == "final") {
    diamond_col <- COL_RE
    pool_label  <- "Meta (RE)"
    stat_line <- sprintf("k=%d  Meta τ²=%.3f  I²=%.0f%%  Qp=%.1e",
                          jrow$k_studies, jrow$meta_tau2, jrow$meta_I2, jrow$meta_Qp)
    eff <- jrow$meta_effect; eff_lo <- jrow$meta_effect_lo; eff_hi <- jrow$meta_effect_hi; eff_p <- jrow$meta_p
  } else {
    diamond_col <- if (kind == "base") COL_BASE else if (kind == "hy1") COL_HY1 else COL_NOT_SELECTED
    pool_label  <- if (kind == "base") "Base Model (RE)" else if (kind == "hy1") "Base + HY1 (RE)" else "Univariate (RE)"
    stat_line   <- sprintf("k=%d  τ²=%.3f  I²=%.0f%%  Qp=%.1e",
                            jrow$k_studies, jrow$tau2, jrow$I2, jrow$p_Q)
    eff <- jrow$effect_RE; eff_lo <- jrow$effect_RE_lo; eff_hi <- jrow$effect_RE_hi; eff_p <- jrow$p_RE
  }

  # Draw stats line
  n_title_lines <- lengths(regmatches(title, gregexpr("\n", title))) + 1
  top_mar <- n_title_lines + 1.4 + 5.5

  sdf <- study_df[study_df$term_key == term_key, ]
  sdf <- sdf[!duplicated(sdf$study), ]
  n_studies <- nrow(sdf)

  if (n_studies == 0) {
    draw_empty_panel(title)
    return(invisible(NULL))
  }

  # Group studies into blocks by cohort disease-stage category (see
  # classify_cohorts())
  sdf$cohort_class <- cohort_class_map[sdf$study]
  sdf$cohort_class[is.na(sdf$cohort_class)] <- "Other"
  sdf <- sdf[order(match(sdf$cohort_class, CATEGORY_ORDER), sdf$study), ]

  # Include the null value (OR/HR = 1) in the range so abline(v = 1) below
  # always falls inside [x_min, x_max] 
  all_lo <- c(sdf$CI_lower, eff_lo, 1)
  all_hi <- c(sdf$CI_upper, eff_hi, 1)
  x_min_data <- min(all_lo, na.rm = TRUE)
  x_max_data <- max(all_hi, na.rm = TRUE)
  x_range    <- max(x_max_data - x_min_data, 0.01)
  x_pad      <- x_range * 0.05
  x_min      <- max(0, x_min_data - x_pad)
  x_max      <- x_max_data + x_pad

  # Bold category-header row before each block of studies
  present_classes <- CATEGORY_ORDER[CATEGORY_ORDER %in% unique(sdf$cohort_class)]
  y_cursor    <- 0
  header_rows <- list()
  study_y     <- integer(n_studies)
  for (cls in present_classes) {
    header_rows[[length(header_rows) + 1]] <- list(y = y_cursor, label = cls)
    y_cursor <- y_cursor + 1
    idx <- which(sdf$cohort_class == cls)
    for (i in idx) {
      study_y[i] <- y_cursor
      y_cursor   <- y_cursor + 1
    }
  }
  n_rows    <- y_cursor
  diamond_y <- n_rows + 0.5

  n_legend <- length(setdiff(CATEGORY_ORDER, "Other"))
  y_lo <- -(n_legend + 2); y_hi <- diamond_y + 1
  line_rows <- par("csi") / ROW_UNIT_IN_SINGLE

  value_strings <- vapply(seq_len(n_studies), function(j)
    format_effect_ci(sdf$effect[j], sdf$CI_lower[j], sdf$CI_upper[j]), character(1))
  re_string <- sprintf("RE: %s", format_effect_ci(eff, eff_lo, eff_hi))
  p_string  <- sprintf("p=%.2e", eff_p)

  study_label_strings <- vapply(seq_len(n_studies), function(j) {
    lbl <- sdf$study[j]
    if (lbl %in% names(COHORT_ANNOTATIONS)) paste0(COHORT_ANNOTATIONS[[lbl]], "  ", lbl) else paste0("  ", lbl)
  }, character(1))
  header_label_strings <- vapply(header_rows, function(h) h$label, character(1))

  legend_label_strings  <- setdiff(CATEGORY_ORDER, "Other")
  legend_label_widths_in <- strwidth(legend_label_strings, units = "inches", cex = 1.0, font = 2)
  names(legend_label_widths_in) <- legend_label_strings

  right_reserve_in <- max(0.5, c(
    strwidth(value_strings, units = "inches", cex = 0.92),
    strwidth(re_string,     units = "inches", cex = 0.88, font = 2),
    strwidth(p_string,      units = "inches", cex = 0.76, font = 2)
  ))
  left_reserve_in <- max(0.5, c(
    strwidth(study_label_strings,  units = "inches", cex = 1.05),
    strwidth(header_label_strings, units = "inches", cex = 1.0, font = 2),
    strwidth(pool_label,           units = "inches", cex = 1.05, font = 2),
    max(legend_label_widths_in) + 0.3
  ))

  line_in   <- par("csi")
  left_mar  <- left_reserve_in  / line_in + 2.5
  right_mar <- right_reserve_in / line_in + 2.5

  max_total_mar_in <- PANEL_WIDTH_IN * 0.6
  total_mar_in <- (left_mar + right_mar) * line_in
  if (total_mar_in > max_total_mar_in) {
    scale     <- max_total_mar_in / total_mar_in
    left_mar  <- left_mar  * scale
    right_mar <- right_mar * scale
  }

  par(mar = c(4, left_mar, top_mar, right_mar))

  plot(NA, xlim = c(x_min, x_max), ylim = c(y_hi, y_lo),
       xlab = spec$effect_axis, ylab = "", yaxt = "n", bty = "l",
       main = title, cex.main = 1.05, font.main = 2)
  mtext(stat_line, side = 3, line = 0.3, cex = 0.72, font = 1, col = "grey25")
  abline(v = 1, col = "grey60", lty = 2, lwd = 0.8)

  usr   <- par("usr")
  x_gap <- (usr[2] - usr[1]) * 0.015

  # Faint divider above every block after the first.
  for (h in header_rows[-1]) abline(h = h$y - 0.5, col = "#eeeeee", lwd = 0.6)

  for (j in seq_len(n_studies)) {
    row   <- sdf[j, ]
    y_j   <- study_y[j]
    col_j <- CATEGORY_COLORS[[row$cohort_class]]

    segments(row$CI_lower, y_j, row$CI_upper, y_j, col = col_j, lwd = 1.8)
    segments(row$CI_lower, y_j - 0.12, row$CI_lower, y_j + 0.12, col = col_j, lwd = 1.8)
    segments(row$CI_upper, y_j - 0.12, row$CI_upper, y_j + 0.12, col = col_j, lwd = 1.8)
    points(row$effect, y_j, pch = 19, col = col_j, cex = 1.0)
    text(usr[2] + x_gap, y_j, value_strings[j],
         adj = 0, cex = 0.92, col = "#222222", xpd = TRUE)
  }

  abline(h = n_rows - 0.2, col = "#cccccc", lwd = 0.8)
  hw <- (eff_hi - eff_lo) / 2
  polygon(x = c(eff - hw, eff, eff + hw, eff),
          y = c(diamond_y, diamond_y + 0.35, diamond_y, diamond_y - 0.35),
          col = diamond_col, border = NA)
  # p-value on its own line below the estimate rather than appended to the
  # same line 
  text(usr[2] + x_gap, diamond_y - line_rows * 0.6, re_string,
       adj = 0, cex = 0.88, font = 2, col = "#32325c", xpd = TRUE)
  text(usr[2] + x_gap, diamond_y + line_rows * 0.9, p_string,
       adj = 0, cex = 0.76, font = 2, col = "#32325c", xpd = TRUE)

  for (h in header_rows) {
    text(x = usr[1] - x_gap, y = h$y, labels = h$label,
         adj = c(1, 0.5), xpd = TRUE, cex = 1.0, font = 2, col = CATEGORY_COLORS[[h$label]])
  }
  for (j in seq_len(n_studies)) {
    text(x = usr[1] - x_gap, y = study_y[j], labels = study_label_strings[j],
         adj = c(1, 0.5), xpd = TRUE, cex = 1.05, font = 1)
  }
  text(x = usr[1] - x_gap, y = diamond_y, labels = pool_label,
       adj = c(1, 0.5), xpd = TRUE, cex = 1.05, font = 2)

  draw_legend(y_lo, legend_label_widths_in)
  invisible(NULL)
}

draw_legend <- function(y_lo, label_widths_in) {
  usr   <- par("usr")
  x_gap <- (usr[2] - usr[1]) * 0.015
  shown_categories <- setdiff(CATEGORY_ORDER, "Other")
  n     <- length(shown_categories)

  in_per_data_x <- par("pin")[1] / diff(usr[1:2])

  ys <- seq(y_lo + 1, by = 1, length.out = n)

  for (i in seq_len(n)) {
    lbl <- shown_categories[i]
    col <- CATEGORY_COLORS[[lbl]]
    text(x = usr[1] - x_gap, y = ys[i], labels = lbl,
         adj = c(1, 0.5), xpd = TRUE, cex = 1.0, font = 2, col = col)
    dot_offset_x <- (label_widths_in[[lbl]] + 0.15) / in_per_data_x
    points(usr[1] - x_gap - dot_offset_x, ys[i],
           pch = 19, cex = 1.0, col = col, xpd = TRUE)
  }
}

# =============================================================================
# MAIN
# =============================================================================

draw_term_png <- function(term_key, lr, cph, cohort_class_map, out_prefix) {
  slug     <- gsub("[^A-Za-z0-9_]", "_", term_key)
  png_path <- file.path(FIGURES_DIR, sprintf("%s_%s_%s.png", out_prefix, slug, TODAY))

  # Size the canvas to this term's own row count rather than a fixed
  # height, so terms with many cohorts don't crowd/overlap and terms with
  # few aren't left sparse. LR and CPH panels share one cell size, so use
  # whichever has more rows.
  n_rows_lr  <- n_term_rows(lr$study$study[lr$study$term_key   == term_key], cohort_class_map)
  n_rows_cph <- n_term_rows(cph$study$study[cph$study$term_key == term_key], cohort_class_map)
  n_rows       <- max(n_rows_lr, n_rows_cph, 1)
  panel_height <- MARGIN_OVERHEAD_IN + (n_rows + LEGEND_RESERVE_UNITS) * ROW_UNIT_IN_SINGLE

  png(png_path, width = PANEL_WIDTH_IN * 2, height = panel_height, units = "in", res = 300, pointsize = 38)
  par(mfrow = c(1, 2), oma = c(2, 0, 2, 0))
  draw_panel(term_key, lr$joint,  lr$study,  MODEL_SPECS$lr,  cohort_class_map)
  draw_panel(term_key, cph$joint, cph$study, MODEL_SPECS$cph, cohort_class_map)
  mtext(format_term(term_key), outer = TRUE, cex = 1.45, font = 2, line = 0.3)
  mtext(COHORT_FOOTNOTE, side = 1, outer = TRUE, cex = 1.05, line = 0.6)
  dev.off()

  cat(sprintf("  Saved: %s\n", basename(png_path)))
  png_path
}

# Main terms (age, sex, disease_duration) as one combined figure.
# Draws from the BASE model (pre-stepwise age + sex + disease duration),
# so these estimates are unadjusted for whatever quadratic/cubic/interaction
# terms stepwise selection added.
draw_main_terms_png <- function(base_lr, base_cph, cohort_class_map) {
  png_path <- file.path(FIGURES_DIR, sprintf("onestage_forest_main_%s.png", TODAY))

  n_rows_all <- unlist(lapply(BASE_TERM_KEYS, function(term_key) {
    c(n_term_rows(base_lr$study$study[base_lr$study$term_key   == term_key], cohort_class_map),
      n_term_rows(base_cph$study$study[base_cph$study$term_key == term_key], cohort_class_map))
  }))
  n_rows       <- max(n_rows_all, 1)
  panel_height <- MARGIN_OVERHEAD_IN + (n_rows + LEGEND_RESERVE_UNITS) * ROW_UNIT_IN_SINGLE

  png(png_path, width = PANEL_WIDTH_IN * 4, height = panel_height * 2 + 1.5, units = "in", res = 300, pointsize = 38)
  par(mfrow = c(2, 4), oma = c(2, 0, 3, 0))
  for (term_key in BASE_TERM_KEYS) {
    draw_panel(term_key, base_lr$joint,  base_lr$study,  MODEL_SPECS$lr,  cohort_class_map, term_label = format_term(term_key), kind = "base")
    draw_panel(term_key, base_cph$joint, base_cph$study, MODEL_SPECS$cph, cohort_class_map, term_label = format_term(term_key), kind = "base")
  }
  plot.new(); plot.new()
  mtext("Main Terms: Age, Sex, Disease Duration (Base Model)", outer = TRUE, cex = 1.5, font = 2, line = 0.8)
  mtext(COHORT_FOOTNOTE, side = 1, outer = TRUE, cex = 1.05, line = 0.6)
  dev.off()

  cat(sprintf("  Saved: %s\n", basename(png_path)))
  png_path
}

# Supplementary figure: hy1_baseline, CPH only
draw_hy1_baseline_png <- function(cph_hy1, cohort_class_map) {
  png_path <- file.path(FIGURES_DIR, sprintf("onestage_forest_hy1_baseline_%s.png", TODAY))

  n_rows       <- max(n_term_rows(cph_hy1$study$study[cph_hy1$study$term_key == "hy1_baseline"], cohort_class_map), 1)
  panel_height <- MARGIN_OVERHEAD_IN + (n_rows + LEGEND_RESERVE_UNITS) * ROW_UNIT_IN_SINGLE

  png(png_path, width = PANEL_WIDTH_IN, height = panel_height, units = "in", res = 300, pointsize = 38)

  par(mfrow = c(1, 1), oma = c(2, 0, 2, 0))
  draw_panel("hy1_baseline", cph_hy1$joint, cph_hy1$study, MODEL_SPECS$cph, cohort_class_map, kind = "hy1")
  mtext(format_term("hy1_baseline"), outer = TRUE, cex = 1.45, font = 2, line = 0.3)
  mtext(COHORT_FOOTNOTE, side = 1, outer = TRUE, cex = 1.05, line = 0.6)
  dev.off()

  cat(sprintf("  Saved: %s\n", basename(png_path)))
  png_path
}

main <- function() {
  cohort_class_map <- classify_cohorts()

  lr  <- load_model_results("lr")
  cph <- load_model_results("cph")

  base_lr  <- load_base_model("lr")
  base_cph <- load_base_model("cph")

  # --- Main figure: base-model terms (age, sex, disease_duration), one
  # combined multi-panel PNG -------------------------------------------
  cat(sprintf("\nMain terms (base model): %s\n", paste(BASE_TERM_KEYS, collapse = ", ")))
  main_file <- draw_main_terms_png(base_lr, base_cph, cohort_class_map)
  cat(sprintf("\nMain comparison figure saved to %s/%s\n", basename(FIGURES_DIR), basename(main_file)))

  # --- Supplementary figures: terms tested during stepwise screening,
  # one PNG per term, unchanged --------------------------------------
  candidate_term_keys <- setdiff(unique(c(lr$joint$term_key, cph$joint$term_key)), BASE_TERM_KEYS)
  cat(sprintf("\nSupplementary (stepwise-tested) terms: %s\n", paste(candidate_term_keys, collapse = ", ")))
  supp_files <- vapply(candidate_term_keys, draw_term_png, character(1),
                        lr = lr, cph = cph, cohort_class_map = cohort_class_map,
                        out_prefix = "onestage_forest_supplementary")
  cat(sprintf("\n%d supplementary plots saved to %s/\n", length(supp_files), basename(FIGURES_DIR)))

  # --- Supplementary figure: hy1_baseline, CPH only. Wrapped in tryCatch
  # since it depends on 05_fit_progression_models.R's optional hy1_baseline
  # step having been run and its output files existing. ---------------------
  tryCatch({
    cph_hy1  <- load_hy1_baseline_model("cph")
    hy1_file <- draw_hy1_baseline_png(cph_hy1, cohort_class_map)
    cat(sprintf("\nHY1-baseline forest plot saved to %s/%s\n", basename(FIGURES_DIR), basename(hy1_file)))
  }, error = function(e) {
    cat(sprintf("\nSkipping hy1_baseline forest plot -- %s\n", e$message))
  })
}

main()
