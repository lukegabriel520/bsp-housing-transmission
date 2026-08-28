#===============================================================================
# 02_model_a_frequentist.R
# MODEL A: Fixed-Effects Panel Interaction Model + Local Projections IRF
#
# SPECIFICATION:
#   Delta ln(RPPI_it) = alpha_i + sum_h [beta_h * (PropType_h x Delta RRP_{t-k})]
#                        + gamma * X_{r,t} + Seasonal_s + Trend_t + eps_it
#
# IDENTIFICATION NOTE:
#   With only N=8 cross-sectional units and one aggregate national policy rate,
#   region/property FE + a single national RRP variable is collinear with time
#   FE (RRP has no cross-sectional variation net of time). We solve this by:
#     (a) interacting Delta RRP with property-type DUMMIES (relative elasticities
#         identified off the panel's cross-sectional heterogeneity), using
#         Single-Detached as the omitted baseline category, and
#     (b) using unit (region x type) FE + quarter-season FE + linear trend,
#         NOT full time FE (which would fully absorb RRP itself).
#   This is standard practice for small-N heterogeneous-panel monetary
#   transmission studies (cf. state-level/sector-level pass-through literature).
#
# INFERENCE: With only 8 clusters, asymptotic cluster-robust SEs are unreliable
#   (Cameron-Gelbach-Miller 2008). We use Wild Cluster Bootstrap with Rademacher
#   weights (999 replications) via fwildclusterboot for valid small-sample
#   inference on every beta_h coefficient.
#===============================================================================

library(tidyverse)
library(fixest)
library(fwildclusterboot)
library(modelsummary)

panel <- read_csv("data/bsp_rppi_panel.csv", show_col_types = FALSE) %>%
  mutate(
    property_type = factor(property_type,
                            levels = c("Single-Detached", "Duplex", "Townhouse", "Condominium")),
    region = factor(region),
    unit_id = factor(unit_id),
    quarter = factor(quarter)
  ) %>%
  filter(!is.na(dln_RPPI), !is.na(d_RRP_lag2))  # drop rows lost to lagging

# -----------------------------------------------------------------------------
# MODEL A1: Baseline interaction model, 2-quarter transmission lag
#   (2 quarters is the modal finding in EM bank-repricing / appraisal-lag
#   literature: banks reprice mortgages with a lag, and RPPI itself reflects
#   loan-closing dates that trail rate decisions by one to two quarters)
# -----------------------------------------------------------------------------
mod_a1 <- feols(
  dln_RPPI ~ property_type : d_RRP_lag2 +
             cpi_national + real_gdp_growth +
             i(quarter) + t | unit_id,
  data = panel,
  cluster = ~unit_id
)

# -----------------------------------------------------------------------------
# MODEL A2: Multi-lag dynamic specification (t-1, t-2, t-4) to trace out timing
# -----------------------------------------------------------------------------
mod_a2 <- feols(
  dln_RPPI ~ property_type : d_RRP_lag1 +
             property_type : d_RRP_lag2 +
             property_type : d_RRP_lag4 +
             cpi_national + real_gdp_growth +
             i(quarter) + t | unit_id,
  data = panel,
  cluster = ~unit_id
)

# -----------------------------------------------------------------------------
# MODEL A3: Lending-rate pass-through channel test
#   (does the RRP -> RPPI effect operate through the bank lending rate?)
# -----------------------------------------------------------------------------
mod_a3 <- feols(
  dln_RPPI ~ property_type : d_lending_rate +
             cpi_national + real_gdp_growth +
             i(quarter) + t | unit_id,
  data = panel,
  cluster = ~unit_id
)

cat("\n============================================================\n")
cat("MODEL A1: Baseline (2-quarter lag) results\n")
cat("============================================================\n")
print(summary(mod_a1))

# -----------------------------------------------------------------------------
# WILD CLUSTER BOOTSTRAP (Rademacher, 999 reps) -- valid small-N inference
#   Run for every property-type interaction coefficient in Model A1
# -----------------------------------------------------------------------------
cat("\n============================================================\n")
cat("WILD CLUSTER BOOTSTRAP (Rademacher weights, 999 reps, N=8 clusters)\n")
cat("============================================================\n")

interaction_coefs <- grep("property_type.*d_RRP_lag2", names(coef(mod_a1)), value = TRUE)

regression_results <- tibble(
  term = interaction_coefs,
  coef = coef(mod_a1)[interaction_coefs],
  se_cluster = se(mod_a1)[interaction_coefs],
  p_value = pvalue(mod_a1)[interaction_coefs]
)
write_csv(regression_results, "data/regression_results.csv")

wcb_results <- map_dfr(interaction_coefs, function(cf) {
  boot <- boottest(
    mod_a1,
    param = cf,
    clustid = "unit_id",
    B = 999,
    type = "rademacher",
    seed = 20260828
  )
  tibble(
    coefficient = cf,
    estimate = boot$point_estimate,
    ci_low = boot$conf_int[1],
    ci_high = boot$conf_int[2],
    p_value_wcb = boot$p_val
  )
})

print(wcb_results)
write_csv(wcb_results, "data/wild_cluster_bootstrap_results.csv")

# -----------------------------------------------------------------------------
# LOCAL PROJECTIONS: 8-quarter horizon IRF of RPPI to a +100bp RRP shock,
# separately by property type (Jorda 2005 LP-IRF, panel version)
# -----------------------------------------------------------------------------
cat("\n============================================================\n")
cat("LOCAL PROJECTIONS: Building 8-quarter horizon IRFs by property type\n")
cat("============================================================\n")

horizons <- 0:8
prop_types <- levels(panel$property_type)

# Build cumulative RPPI log-level (h-step-ahead) as LHS for each horizon
panel_lp <- panel %>%
  arrange(unit_id, t) %>%
  group_by(unit_id) %>%
  mutate(ln_RPPI = log(RPPI)) %>%
  ungroup()

lp_results <- expand_grid(h = horizons, property_type = prop_types) %>%
  mutate(beta = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_)

for (h in horizons) {
  panel_h <- panel_lp %>%
    group_by(unit_id) %>%
    arrange(t, .by_group = TRUE) %>%
    mutate(y_h = lead(ln_RPPI, h) - ln_RPPI) %>%  # cumulative log change over h quarters
    ungroup() %>%
    filter(!is.na(y_h), !is.na(d_RRP_lag1))

  mod_h <- tryCatch(
    feols(
      y_h ~ property_type : d_RRP_lag1 + cpi_national + real_gdp_growth | unit_id,
      data = panel_h,
      cluster = ~unit_id
    ),
    error = function(e) NULL
  )

  if (!is.null(mod_h)) {
    cf <- coef(mod_h)
    se <- se(mod_h)
    for (pt in prop_types) {
      nm <- grep(paste0("property_type", gsub("([()\\-])", "\\\\\\1", pt), ":d_RRP_lag1"),
                 names(cf), value = TRUE)
      if (length(nm) == 1) {
        idx <- which(lp_results$h == h & lp_results$property_type == pt)
        lp_results$beta[idx] <- cf[nm]
        lp_results$se[idx] <- se[nm]
        lp_results$ci_low[idx] <- cf[nm] - 1.96 * se[nm]
        lp_results$ci_high[idx] <- cf[nm] + 1.96 * se[nm]
      }
    }
  }
}

# Convert per-100bp: d_RRP is already in percentage points (pp), so a +1.00
# coefficient on d_RRP_lag1 already represents effect of +100bp. Scale IRF
# to cumulative % price response to a +100bp (1.00 pp) hike:
lp_results <- lp_results %>%
  mutate(
    irf_pct = beta * 100,        # convert log-change coefficient to approx %
    irf_ci_low = ci_low * 100,
    irf_ci_high = ci_high * 100
  )

write_csv(lp_results, "data/local_projections_irf.csv")
cat("Local projections IRF saved to data/local_projections_irf.csv\n")

# -----------------------------------------------------------------------------
# PUBLICATION-GRADE REGRESSION TABLE (modelsummary)
# -----------------------------------------------------------------------------
models_list <- list(
  "2Q Lag"        = mod_a1,
  "Multi-Lag"     = mod_a2,
  "Lending-Rate"  = mod_a3
)

modelsummary(
  models_list,
  stars = c('*' = .1, '**' = .05, '***' = .01),
  gof_omit = "IC|Log|Adj|Within",
  title = "Table 1. BSP Policy Rate Transmission to Residential Property Prices (Fixed-Effects Panel)",
  notes = "Standard errors clustered by unit (region x property type), N=8 clusters. Wild cluster bootstrap p-values reported separately due to small cluster count (see wild_cluster_bootstrap_results.csv).",
  output = "outputs/table1_regression_results.docx"
)

modelsummary(
  models_list,
  stars = c('*' = .1, '**' = .05, '***' = .01),
  gof_omit = "IC|Log|Adj|Within",
  output = "outputs/table1_regression_results.txt"
)

saveRDS(list(mod_a1 = mod_a1, mod_a2 = mod_a2, mod_a3 = mod_a3),
        "data/model_a_objects.rds")

cat("\nModel A complete. Regression table saved to outputs/table1_regression_results.docx\n")
