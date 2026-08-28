#===============================================================================
# 04_diagnostics.R
# MODULE 3: Panel Diagnostics & Specification Tests
#
# 1. Panel unit-root tests: Levin-Lin-Chu (common unit root) and
#    Im-Pesaran-Shin (individual unit roots) on RPPI log-levels and
#    RRP -- confirms whether we need to work in differences (we do:
#    dln_RPPI is I(0) by construction as a log-difference).
# 2. Wooldridge test for serial correlation in panel data.
# 3. Breusch-Pagan / Studentized BP test for heteroskedasticity.
# 4. Pass-through completeness validation: is RRP -> lending rate
#    transmission complete (coefficient = 1) or partial in the long run?
#===============================================================================

library(tidyverse)
library(plm)
library(lmtest)

panel_raw <- read_csv("data/bsp_rppi_panel.csv", show_col_types = FALSE)

pdata <- pdata.frame(panel_raw, index = c("unit_id", "t"))

cat("============================================================\n")
cat("1. PANEL UNIT ROOT TESTS\n")
cat("============================================================\n")

# --- On RPPI log-levels (expect: unit root present, i.e. I(1)) ---
cat("\n--- ln(RPPI) level series ---\n")
pdata$ln_RPPI <- log(pdata$RPPI)

llc_level <- tryCatch(
  purtest(ln_RPPI ~ 1, data = pdata, index = c("unit_id", "t"),
          test = "levinlin", lags = "SIC"),
  error = function(e) { cat("LLC error:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(llc_level)) {
  cat("Levin-Lin-Chu (common unit root), ln(RPPI):\n")
  print(llc_level)
}

ips_level <- tryCatch(
  purtest(ln_RPPI ~ 1, data = pdata, index = c("unit_id", "t"),
          test = "ips", lags = "SIC"),
  error = function(e) { cat("IPS error:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(ips_level)) {
  cat("\nIm-Pesaran-Shin (individual unit roots), ln(RPPI):\n")
  print(ips_level)
}

# --- On dln_RPPI (first-differenced; expect: stationary, reject unit root) ---
cat("\n--- Delta ln(RPPI), first-differenced series ---\n")
pdata_d <- pdata[!is.na(pdata$dln_RPPI), ]

llc_diff <- tryCatch(
  purtest(dln_RPPI ~ 1, data = pdata_d, index = c("unit_id", "t"),
          test = "levinlin", lags = "SIC"),
  error = function(e) { cat("LLC error:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(llc_diff)) {
  cat("Levin-Lin-Chu, Delta ln(RPPI):\n")
  print(llc_diff)
}

ips_diff <- tryCatch(
  purtest(dln_RPPI ~ 1, data = pdata_d, index = c("unit_id", "t"),
          test = "ips", lags = "SIC"),
  error = function(e) { cat("IPS error:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(ips_diff)) {
  cat("\nIm-Pesaran-Shin, Delta ln(RPPI):\n")
  print(ips_diff)
}

cat("\nINTERPRETATION: If level series fail to reject the unit-root null\n")
cat("(non-stationary, I(1)) while the differenced series rejects it\n")
cat("(stationary, I(0)), this CONFIRMS our modeling choice of working with\n")
cat("Delta ln(RPPI) as the dependent variable throughout Modules 2A/2B,\n")
cat("avoiding spurious regression in levels.\n")

# -----------------------------------------------------------------------------
cat("\n============================================================\n")
cat("2. SERIAL CORRELATION: Wooldridge Test for Panel Data\n")
cat("============================================================\n")

plm_fe <- plm(
  dln_RPPI ~ d_RRP_lag2 + cpi_national + real_gdp_growth,
  data = pdata, model = "within"
)

wooldridge_test <- tryCatch(
  pwartest(plm_fe),
  error = function(e) { cat("Wooldridge test error:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(wooldridge_test)) print(wooldridge_test)

cat("\nINTERPRETATION: p < 0.05 rejects the null of no serial correlation,\n")
cat("indicating within-unit residuals are autocorrelated -- justifies our\n")
cat("use of cluster-robust (and wild-bootstrap) SEs in Model A rather than\n")
cat("classical iid SEs.\n")

# -----------------------------------------------------------------------------
cat("\n============================================================\n")
cat("3. HETEROSKEDASTICITY: Breusch-Pagan / Studentized BP Test\n")
cat("============================================================\n")

bp_test <- tryCatch(
  bptest(plm_fe),
  error = function(e) { cat("BP test error:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(bp_test)) print(bp_test)

cat("\nINTERPRETATION: p < 0.05 rejects homoskedasticity -- further supports\n")
cat("cluster-robust / heteroskedasticity-robust inference throughout.\n")

# -----------------------------------------------------------------------------
cat("\n============================================================\n")
cat("4. PASS-THROUGH COMPLETENESS: BSP RRP -> Commercial Lending Rate\n")
cat("============================================================\n")

pt_data <- panel_raw %>%
  distinct(t, RRP, lending_rate, d_RRP, d_lending_rate) %>%
  arrange(t) %>%
  filter(!is.na(d_RRP), !is.na(d_lending_rate))

pt_model <- lm(d_lending_rate ~ d_RRP, data = pt_data)
cat("\nContemporaneous pass-through regression: Delta(lending_rate) ~ Delta(RRP)\n")
print(summary(pt_model))

pt_coef <- coef(pt_model)["d_RRP"]
pt_se <- summary(pt_model)$coefficients["d_RRP", "Std. Error"]

cat("\nContemporaneous pass-through coefficient:", round(pt_coef, 3), "\n")
cat("(1.0 = complete/full pass-through, <1.0 = partial pass-through)\n")

# Test H0: pass-through = 1 (complete pass-through)
t_stat_complete <- (pt_coef - 1) / pt_se
p_val_complete <- 2 * pt(-abs(t_stat_complete), df = pt_model$df.residual)
cat("\nTest of H0: pass-through coefficient = 1 (complete pass-through)\n")
cat("  t-statistic:", round(t_stat_complete, 3), "\n")
cat("  p-value:", round(p_val_complete, 4), "\n")
if (p_val_complete < 0.05) {
  cat("  => REJECT H0: pass-through is statistically INCOMPLETE in the same quarter.\n")
} else {
  cat("  => FAIL TO REJECT H0: consistent with complete pass-through.\n")
}

# Long-run pass-through via cumulative lending rate adjustment (partial adj model)
cat("\n--- Long-run pass-through (cumulative over 4 quarters) ---\n")
pt_data_lr <- panel_raw %>%
  distinct(t, RRP, lending_rate) %>%
  arrange(t) %>%
  mutate(
    RRP_chg_4q = RRP - lag(RRP, 4),
    lend_chg_4q = lending_rate - lag(lending_rate, 4)
  ) %>%
  filter(!is.na(RRP_chg_4q), !is.na(lend_chg_4q))

pt_model_lr <- lm(lend_chg_4q ~ RRP_chg_4q, data = pt_data_lr)
print(summary(pt_model_lr))
cat("\nLong-run (4-quarter) pass-through coefficient:",
    round(coef(pt_model_lr)["RRP_chg_4q"], 3), "\n")

# -----------------------------------------------------------------------------
# SAVE DIAGNOSTIC SUMMARY
# -----------------------------------------------------------------------------
diag_summary <- tibble(
  test = c("Wooldridge serial correlation", "Breusch-Pagan heteroskedasticity",
           "Contemporaneous pass-through", "Long-run (4Q) pass-through"),
  statistic = c(
    if (!is.null(wooldridge_test)) wooldridge_test$statistic else NA,
    if (!is.null(bp_test)) bp_test$statistic else NA,
    pt_coef,
    coef(pt_model_lr)["RRP_chg_4q"]
  ),
  p_value = c(
    if (!is.null(wooldridge_test)) wooldridge_test$p.value else NA,
    if (!is.null(bp_test)) bp_test$p.value else NA,
    summary(pt_model)$coefficients["d_RRP", "Pr(>|t|)"],
    summary(pt_model_lr)$coefficients["RRP_chg_4q", "Pr(>|t|)"]
  )
)
write_csv(diag_summary, "data/diagnostic_tests.csv")
cat("\nDiagnostic summary saved to data/diagnostic_tests.csv\n")
