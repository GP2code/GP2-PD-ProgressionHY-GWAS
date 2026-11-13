library(tidyverse)
library(survival)
library(survminer)
library(broom)
library(yaml)

# ── Load config ───────────────────────────────────────────────────────────────
cfg <- read_yaml("config.yaml")

# ============================================================
# LOGGING SETUP
# ============================================================

dir.create("logs",    showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
log_file <- cfg$analysis_log
cat(paste("=== HY3 Analysis Log ===\n"), file = log_file)
cat(paste("Analysis started:", Sys.time(), "\n\n"), file = log_file, append = TRUE)

log_msg <- function(msg, console = TRUE) {
  if (console) cat(msg)
  cat(msg, file = log_file, append = TRUE)
}

# Helper to log print() output (tables, model summaries, tibbles)
log_print <- function(x, console = TRUE) {
  out <- paste(capture.output(
    if (inherits(x, c("tbl", "tbl_df", "data.frame"))) {
      print(x, n = Inf, width = Inf)
    } else {
      print(x)
    }
  ), collapse = "\n")
  out <- paste0(out, "\n")
  if (console) cat(out)
  cat(out, file = log_file, append = TRUE)
}

# ============================================================
# 0. LOAD DATA
# ============================================================
on_df  <- read_csv(cfg$on_clinic_path)
off_df <- read_csv(cfg$off_clinic_path)

# Combined for distribution comparison
both_df <- bind_rows(on_df, off_df)

# ============================================================
# 1. HY OFF vs ON CROSSTAB (TWO-ASSESSMENT VISITS)
# ============================================================
raw_df <- read_csv(cfg$raw_hy_path)

# Keep only visits with exactly 2 assessments (any state)
two_assessment <- raw_df |>
  group_by(id, visit_name) |>
  summarise(n_assessments = n(), .groups = "drop") |>
  filter(n_assessments == 2) |>
  select(id, visit_name)

df_kept <- raw_df |> inner_join(two_assessment, by = c("id", "visit_name"))

log_msg(sprintf("Two-assessment visits: %d\n", nrow(two_assessment)))
log_msg(sprintf("Unique IDs:            %d\n", df_kept |> distinct(id) |> nrow()))

# Split into ON and OFF
on_assess  <- df_kept |> filter(state == "ON")  |> select(id, visit_name, hy_on  = hy)
off_assess <- df_kept |> filter(state == "OFF") |> select(id, visit_name, hy_off = hy)

paired <- off_assess |>
  inner_join(on_assess, by = c("id", "visit_name")) |>
  filter(!is.na(hy_off), !is.na(hy_on))

# Crosstab HY OFF vs HY ON
log_msg("\n--- HY OFF vs ON crosstab (two-assessment visits) ---\n")
crosstab <- table("HY (OFF)" = paired$hy_off, "HY (ON)" = paired$hy_on)
log_print(addmargins(crosstab))

# Count two-assessment visits with Unknown in either assessment
unknown_two <- raw_df |>
  inner_join(two_assessment, by = c("id", "visit_name")) |>
  group_by(id, visit_name) |>
  summarise(states = list(unique(state)), .groups = "drop") |>
  filter(map_lgl(states, ~ "UNKNOWN" %in% .x & any(c("ON", "OFF") %in% .x)))

log_msg(sprintf("\nTwo-assessment visits with Unknown in either assessment: %d\n",
            nrow(unknown_two)))
log_msg(sprintf("Unique IDs with such visits: %d\n",
            unknown_two |> distinct(id) |> nrow()))

# ============================================================
# 2. SAMPLE PIVOT TABLE (3 IDs WITH >= 1 ON+OFF VISIT)
# ============================================================

# IDs with at least one visit assessed in both ON and OFF
multi_state_ids <- raw_df |>
  group_by(id, visit_name) |>
  summarise(states = list(unique(state)), .groups = "drop") |>
  filter(map_lgl(states, ~ "ON" %in% .x & "OFF" %in% .x)) |>
  distinct(id) |>
  pull(id)

log_msg(sprintf("Participants with >= 1 visit assessed in both ON and OFF: %d\n",
            length(multi_state_ids)))

# Sample 3 IDs reproducibly
set.seed(42)
sample_ids <- sample(multi_state_ids, size = min(3, length(multi_state_ids)), replace = FALSE)

# All visits for sampled IDs, pivot one row per id+visit_name
visit_order <- c("BL", paste0("V", sprintf("%02d", 1:22)))

pivot_table <- raw_df |>
  filter(id %in% sample_ids) |>
  group_by(id, visit_name, visit_month, med_for_pd) |>
  reframe(
    hy_ON          = hy[state == "ON"][1],
    hy_OFF         = hy[state == "OFF"][1],
    hy_Unknown     = hy[state == "Unknown"][1],
    updrs3_ON      = updrs3[state == "ON"][1],
    updrs3_OFF     = updrs3[state == "OFF"][1],
    updrs3_Unknown = updrs3[state == "Unknown"][1]
  ) |>
  distinct() |>
  arrange(id, visit_month) |>
  select(id, visit_month, med_for_pd,
         hy_ON, hy_OFF, hy_Unknown,
         updrs3_ON, updrs3_OFF, updrs3_Unknown)

log_msg("\n--- Sample pivot table (3 IDs, all visits) ---\n")
log_print(pivot_table)

# ============================================================
# 4. BASELINE STATE COUNTS DIAGNOSTIC
#    Note: most BL visits are Unknown state, so logistic
#    regression results will be near-identical for both clinics
# ============================================================

log_msg("\n--- Baseline state counts (visit_month == 0) ---\n")
on_df |> filter(visit_month == 0) |> count(clinic)
off_df |> filter(visit_month == 0) |> count(clinic)

# ============================================================
# 2. HY DISTRIBUTION COMPARISON
# ============================================================

log_msg("\n--- HY Distribution: On Clinic vs Off Clinic ---\n")

dist_table <- both_df |>
  filter(!is.na(hy)) |>
  count(clinic, hy) |>
  group_by(clinic) |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  pivot_wider(names_from = clinic, values_from = c(n, pct))

log_print(dist_table)

# Chi-square test
chi_mat <- both_df |>
  filter(!is.na(hy)) |>
  count(clinic, hy) |>
  pivot_wider(names_from = hy, values_from = n, values_fill = 0) |>
  column_to_rownames("clinic") |>
  as.matrix()

chi_result <- chisq.test(chi_mat)
log_msg(sprintf("\nChi-square: X²=%.2f, df=%d, p=%.4f\n",
            chi_result$statistic, chi_result$parameter, chi_result$p.value))

# ============================================================
# 3. LOGISTIC REGRESSION AT BASELINE
#    Note: results near-identical across clinics as most of
#    baseline visits have Unknown state.
# ============================================================

log_msg("\n--- Logistic Regression at BL: HY3 ~ age + sex + disease_duration ---\n")

logit_results <- list()
for (clinic in c("On Clinic", "Off Clinic")) {
  bl <- both_df |>
    filter(clinic == !!clinic, visit_month == 0) |>
    filter(!is.na(hy3), !is.na(age), !is.na(sex), !is.na(disease_duration)) |>
    mutate(sex = factor(sex, levels = c("Female", "Male")))

  log_msg(sprintf("\n%s — N=%d, HY3 prevalence=%.1f%%\n",
              clinic, nrow(bl), 100 * mean(bl$hy3)))

  model <- glm(hy3 ~ age + sex + disease_duration,
               data = bl, family = binomial())

  res <- tidy(model, conf.int = TRUE)
  log_print(res |> select(term, estimate, std.error, statistic, p.value, conf.low, conf.high))
  logit_results[[clinic]] <- res
}

# ============================================================
# 4. SURVIVAL ANALYSIS: TIME TO HY >= 3
#    - Exclude IDs with HY >= 3 at baseline
#    - Event: first visit with HY >= 3 after baseline
#    - Censored: last visit with HY < 3
# ============================================================

log_msg("\n--- Survival Analysis: Time to HY >= 3 ---\n")

build_survival_data <- function(df) {
  # Exclude IDs with HY >= 3 at baseline
  bl_hy <- df |>
    filter(visit_month == 0) |>
    group_by(id) |>
    summarise(bl_hy = first(hy), .groups = "drop")

  eligible_ids <- bl_hy |> filter(bl_hy < 3) |> pull(id)

  baseline_covars <- df |>
    filter(id %in% eligible_ids, visit_month == 0) |>
    group_by(id) |>
    summarise(age = first(age), sex = first(sex),
              disease_duration = first(disease_duration), .groups = "drop")

  # Find event/censoring time per ID
  records <- df |>
    filter(id %in% eligible_ids, visit_month > 0) |>
    arrange(id, visit_month) |>
    group_by(id) |>
    summarise(
      time = case_when(
        any(hy >= 3, na.rm = TRUE) ~ first(visit_month[hy >= 3 & !is.na(hy)]),
        TRUE                        ~ last(visit_month[hy < 3  & !is.na(hy)])
      ),
      event = as.integer(any(hy >= 3, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    filter(!is.na(time), time > 0)

  records |>
    inner_join(baseline_covars, by = "id") |>
    mutate(sex = factor(sex, levels = c("Female", "Male")),
           sex_male = as.integer(sex == "Male")) |>
    drop_na(time, event, age, sex_male, disease_duration)
}

surv_data <- list(
  "On Clinic"  = build_survival_data(on_df),
  "Off Clinic" = build_survival_data(off_df)
)

cox_results <- list()
for (clinic in names(surv_data)) {
  sd <- surv_data[[clinic]]
  log_msg(sprintf("\n%s — N=%d, Events=%d\n", clinic, nrow(sd), sum(sd$event)))
  cox <- coxph(Surv(time, event) ~ age + sex_male + disease_duration, data = sd)
  cox_results[[clinic]] <- tidy(cox, conf.int = TRUE, exponentiate = FALSE)
  log_print(cox_results[[clinic]] |> select(term, estimate, std.error, statistic, p.value,
                                         conf.low, conf.high))
}

# ============================================================
# 5. PLOTS
# ============================================================

colors <- c("On Clinic" = "#2196F3", "Off Clinic" = "#E91E63")

# --- Panel A: KM curves (both clinics on same plot) ---
combined_surv <- bind_rows(
  surv_data[["On Clinic"]]  |> mutate(clinic = "On Clinic"),
  surv_data[["Off Clinic"]] |> mutate(clinic = "Off Clinic")
)

combined_surv <- combined_surv |>
  mutate(clinic = factor(clinic, levels = c("On Clinic", "Off Clinic")))

km_fit <- survfit(Surv(time, event) ~ clinic, data = combined_surv)

# Log-rank test -- logged explicitly in addition to being displayed on the
# plot itself via ggsurvplot's pval argument below.
logrank <- survdiff(Surv(time, event) ~ clinic, data = combined_surv)
logrank_p <- pchisq(logrank$chisq, df = length(logrank$n) - 1, lower.tail = FALSE)
log_msg(sprintf("\nLog-rank test (On Clinic vs Off Clinic): chisq=%.3f, df=%d, p=%.4f\n",
            logrank$chisq, length(logrank$n) - 1, logrank_p))

p_km <- ggsurvplot(
  km_fit,
  data        = combined_surv,
  palette     = c("#2196F3", "#E91E63"),  # On Clinic, Off Clinic
  conf.int    = TRUE,
  legend.labs = c("On Clinic", "Off Clinic"),
  pval        = TRUE,
  pval.method = TRUE,
  pval.size   = 4,
  xlab        = "Months from baseline",
  ylab        = "Survival probability",
  title       = "Time to HY >= 3: On Clinic vs Off Clinic",
  ggtheme     = theme_classic(base_size = 12),
  risk.table  = TRUE,
  risk.table.height = 0.25
)

# --- Panel B: Logistic regression forest plot ---
logit_forest_df <- bind_rows(
  logit_results[["On Clinic"]]  |> mutate(clinic = "On Clinic"),
  logit_results[["Off Clinic"]] |> mutate(clinic = "Off Clinic")
) |>
  filter(term != "(Intercept)") |>
  mutate(
    term = recode(term,
      "age"             = "Age",
      "sexMale"         = "Sex (Male vs Female)",
      "disease_duration"= "Disease Duration"
    ),
    term = factor(term, levels = c("Disease Duration", "Sex (Male vs Female)", "Age"))
  )

p_logit <- ggplot(logit_forest_df,
                  aes(x = estimate, y = term, color = clinic,
                      xmin = conf.low, xmax = conf.high)) +
  geom_pointrange(position = position_dodge(width = 0.4), size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = colors) +
  labs(x = "Log Odds (95% CI)", y = NULL, color = NULL,
       title = "Logistic Regression at BL",
       subtitle = "HY \u2265 3 ~ age + sex + disease_duration") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

# --- Panel C: Cox PH forest plot ---
cox_forest_df <- bind_rows(
  cox_results[["On Clinic"]]  |> mutate(clinic = "On Clinic"),
  cox_results[["Off Clinic"]] |> mutate(clinic = "Off Clinic")
) |>
  mutate(
    term = recode(term,
      "age"             = "Age",
      "sex_male"        = "Sex (Male vs Female)",
      "disease_duration"= "Disease Duration"
    ),
    term = factor(term, levels = c("Disease Duration", "Sex (Male vs Female)", "Age"))
  )

p_cox <- ggplot(cox_forest_df,
                aes(x = estimate, y = term, color = clinic,
                    xmin = conf.low, xmax = conf.high)) +
  geom_pointrange(position = position_dodge(width = 0.4), size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = colors) +
  labs(x = "Log Hazard Ratio (95% CI)", y = NULL, color = NULL,
       title = "Cox PH: Time to HY \u2265 3",
       subtitle = "Hazard ~ age + sex + disease_duration") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

# ============================================================
# 6. HY TRAJECTORY PLOT
#    Mean HY ± SEM per visit_month, by state (ON / OFF / Unknown)
#    Uses raw df passed from wrangling — read separately if needed
# ============================================================

# Read the raw (pre-clinic-split) data for trajectory
# This needs all three states, so we use the original df saved from Python
# raw_df already loaded in Section 1

traj_colors <- c("ON" = "#2196F3", "OFF" = "#E91E63", "UNKNOWN" = "#9C27B0")
traj_labels <- c("ON" = "ON", "OFF" = "OFF", "UNKNOWN" = "Unknown")

# Summarise HY per visit_month + state
traj_df <- raw_df |>
  filter(!is.na(hy), !is.na(visit_month), state %in% c("ON", "OFF", "UNKNOWN")) |>
  group_by(visit_month, state) |>
  summarise(
    mean_hy = mean(hy, na.rm = TRUE),
    sem_hy  = sd(hy, na.rm = TRUE) / sqrt(n()),
    n       = n(),
    .groups = "drop"
  ) |>
  filter(n >= 5) |>
  arrange(state, visit_month)

x_breaks <- seq(0, max(traj_df$visit_month, na.rm = TRUE), by = 6)

p_traj <- ggplot(traj_df, aes(x = visit_month, y = mean_hy,
                               color = state, linetype = state,
                               shape = state, group = state)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  geom_ribbon(aes(ymin = mean_hy - sem_hy, ymax = mean_hy + sem_hy, fill = state),
              alpha = 0.12, color = NA) +
  scale_color_manual(values = traj_colors, labels = traj_labels) +
  scale_fill_manual(values  = traj_colors, labels = traj_labels) +
  scale_linetype_manual(values = c("ON" = "solid", "OFF" = "dashed", "UNKNOWN" = "dotted"),
                        labels = traj_labels) +
  scale_shape_manual(values = c("ON" = 16, "OFF" = 15, "UNKNOWN" = 17),
                     labels = traj_labels) +
  scale_x_continuous(breaks = x_breaks) +
  labs(
    x = "Months from baseline", y = "Mean H&Y stage",
    title = "Mean H&Y Stage Trajectory by Medication State",
    color = "State", fill = "State", linetype = "State", shape = "State"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")

# Sample size table
n_table <- traj_df |>
  select(visit_month, state, n) |>
  pivot_wider(names_from = state, values_from = n, values_fill = 0) |>
  arrange(visit_month)

log_msg("\n--- Sample sizes per visit and state ---\n")
log_print(n_table)

# --- Save plots ---
# KM saved via png device (ggsurvplot objects don't work with ggsave)
png(cfg$km_plot, width = 10, height = 7, units = "in", res = 150)
print(p_km)
dev.off()
 
# Forest plots side by side
forest_combined <- cowplot::plot_grid(p_logit, p_cox, ncol = 2, labels = c("B", "C"))
ggsave(cfg$forest_plots, plot = forest_combined, width = 12, height = 5, dpi = 150)
 
# Trajectory plot
ggsave(cfg$trajectory_plot, plot = p_traj, width = 12, height = 4, dpi = 150)
 
log_msg(paste0("\nPlots saved: ", cfg$km_plot, ", ", cfg$forest_plots, ", ", cfg$trajectory_plot, "\n"))