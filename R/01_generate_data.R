#===============================================================================
# 01_generate_data.R
# BSP Monetary Policy -> Philippine Residential Real Estate Transmission Study
#
# PURPOSE: Build a reproducible, empirically-calibrated SYNTHETIC panel dataset
#          that matches the statistical properties (means, volatility, turning
#          points) of actual BSP/PSA published series, because raw microdata
#          from BSP's loan-level RPPI methodology is not public at the granular
#          quarterly x region x property-type cell level.
#
# PANEL STRUCTURE: 8 quarters/type of NCR|AONCR x 4 property types
#                   Q1 2018 - Q2 2026 (34 quarters) => N = 8, T = 34, obs = 272
#
# CALIBRATION ANCHORS (real, sourced Aug 2026):
#   - RRP path: 3.00% (2018 start) -> 4.75% (late 2018) -> 4.00% (2019 cuts)
#               -> 2.00% (2020-21 pandemic floor) -> 6.50% (Oct 2023 peak)
#               -> steady cuts through 2024-25 -> 4.25% (Feb 19, 2026)
#               -> surprise hikes to 4.50% (Apr 2026) and 4.75% (Jun 2026)
#                 on an oil/fertilizer-driven inflation shock
#   - RPPI: NCR +13.9% YoY (Q1 2025), nationwide decel to +1.6% YoY (Q4 2025,
#           6-yr low), rebound to +4.5% YoY / +5.6% QoQ (Q1 2026) led by condos
#   - Condo prices materially more volatile / rate-sensitive than horizontal
#     housing (single-detached), consistent with investor-driven demand
#     and higher loan-to-value / adjustable financing exposure
#===============================================================================

library(tidyverse)
library(lubridate)

set.seed(20260828)  # reproducibility seed = run date, YYYYMMDD

# -----------------------------------------------------------------------------
# 1. QUARTERLY TIME INDEX
# -----------------------------------------------------------------------------
quarters_seq <- seq(as.Date("2018-01-01"), as.Date("2026-04-01"), by = "3 months")
n_q <- length(quarters_seq)  # 34 quarters

time_df <- tibble(
  date    = quarters_seq,
  year    = year(date),
  quarter = quarter(date),
  t       = row_number(),
  yq_lab  = paste0(year, "Q", quarter)
)

# -----------------------------------------------------------------------------
# 2. BSP TARGET RRP RATE PATH (real, piecewise-anchored + light AR noise)
#    Hand-coded from actual Monetary Board decisions (see source comments)
# -----------------------------------------------------------------------------
rrp_anchors <- tribble(
  ~date,          ~rrp,
  "2018-01-01",   3.00,   # start of hiking cycle vs 2018 inflation surge
  "2018-04-01",   3.25,
  "2018-07-01",   4.00,
  "2018-10-01",   4.75,   # 2018 peak - inflation surge response
  "2019-01-01",   4.75,
  "2019-04-01",   4.50,
  "2019-07-01",   4.25,
  "2019-10-01",   4.00,   # 2019 easing
  "2020-01-01",   3.75,
  "2020-04-01",   2.75,   # COVID emergency cuts begin
  "2020-07-01",   2.25,
  "2020-10-01",   2.00,   # pandemic floor
  "2021-01-01",   2.00,
  "2021-04-01",   2.00,
  "2021-07-01",   2.00,
  "2021-10-01",   2.00,
  "2022-01-01",   2.00,
  "2022-04-01",   2.50,   # tightening cycle begins (May 2022)
  "2022-07-01",   3.75,
  "2022-10-01",   5.00,   # +75bp Nov 2022 outsized hike
  "2023-01-01",   5.75,
  "2023-04-01",   6.25,
  "2023-07-01",   6.25,
  "2023-10-01",   6.50,   # cycle peak (Oct 2023)
  "2024-01-01",   6.50,
  "2024-04-01",   6.50,
  "2024-07-01",   6.25,   # easing begins Aug 2024
  "2024-10-01",   6.00,
  "2025-01-01",   5.75,
  "2025-04-01",   5.50,
  "2025-07-01",   5.25,
  "2025-10-01",   5.00,
  "2026-01-01",   4.50,   # includes Feb 19 2026 cut to 4.25% (Q1 avg)
  "2026-04-01",   4.625   # blends Apr hike to 4.50% and Jun hike to 4.75%
) %>%
  mutate(date = as.Date(date))

rrp_df <- time_df %>%
  left_join(rrp_anchors, by = "date") %>%
  mutate(
    RRP = rrp,
    d_RRP = RRP - lag(RRP)  # quarter-on-quarter change in policy rate (pp)
  ) %>%
  select(-rrp)

# -----------------------------------------------------------------------------
# 3. MACRO CONTROLS: Regional CPI inflation & Real GDP growth
#    Calibrated to broad PSA/PSA-BSP published ranges (illustrative AR(1))
# -----------------------------------------------------------------------------
simulate_ar1 <- function(n, mean, sd, rho, start = mean) {
  x <- numeric(n)
  x[1] <- start
  for (i in 2:n) {
    x[i] <- mean + rho * (x[i-1] - mean) + rnorm(1, 0, sd)
  }
  x
}

macro_df <- time_df %>%
  mutate(
    # headline CPI inflation, YoY % -- captures 2018 surge (~5-6%), 2020-21 low
    # (~2-3%), 2022-23 spike (peak ~8.7% in Jan 2023), 2024-25 moderation,
    # 2026 oil/fertilizer-driven re-acceleration toward BSP's revised ~6.3-6.4%
    # full-year forecast
    cpi_national = case_when(
      year == 2018                     ~ simulate_ar1(n(), 5.2, 0.3, 0.6, 4.0)[t],
      year %in% c(2019)                ~ simulate_ar1(n(), 2.6, 0.3, 0.6, 3.8)[t],
      year %in% c(2020, 2021)          ~ simulate_ar1(n(), 2.6, 0.4, 0.6, 2.2)[t],
      year == 2022                     ~ simulate_ar1(n(), 5.8, 0.5, 0.6, 3.0)[t],
      year == 2023                     ~ simulate_ar1(n(), 6.6, 0.7, 0.6, 8.0)[t],
      year == 2024                     ~ simulate_ar1(n(), 3.6, 0.4, 0.6, 3.9)[t],
      year == 2025                     ~ simulate_ar1(n(), 2.0, 0.3, 0.6, 1.8)[t],
      year == 2026                     ~ simulate_ar1(n(), 6.0, 0.5, 0.6, 4.5)[t],
      TRUE                             ~ 3.0
    ),
    real_gdp_growth = case_when(
      year %in% c(2020)                ~ simulate_ar1(n(), -8.0, 2.0, 0.5, -0.7)[t],
      year %in% c(2021)                ~ simulate_ar1(n(), 6.5, 1.5, 0.5, 3.5)[t],
      year %in% c(2022)                ~ simulate_ar1(n(), 7.0, 1.0, 0.5, 8.2)[t],
      year %in% c(2023, 2024)          ~ simulate_ar1(n(), 5.7, 0.6, 0.5, 5.5)[t],
      year %in% c(2025, 2026)          ~ simulate_ar1(n(), 5.2, 0.6, 0.5, 5.0)[t],
      TRUE                              ~ simulate_ar1(n(), 6.3, 0.5, 0.5, 6.3)[t]
    )
  ) %>%
  select(date, t, cpi_national, real_gdp_growth)

# -----------------------------------------------------------------------------
# 4. COMMERCIAL BANK LENDING RATE PASS-THROUGH
#    Universal/commercial bank lending rate = markup over RRP with partial,
#    lagged pass-through (banks reprice with delay -> incomplete transmission
#    in short run, converging toward long-run pass-through ~0.7-0.9)
# -----------------------------------------------------------------------------
pt_df <- rrp_df %>%
  arrange(t) %>%
  mutate(
    # Long-run spread of ~3.1pp over RRP for housing loans (illustrative,
    # consistent with BSP published UCB lending rate series behavior)
    base_spread = 3.1,
    # Partial adjustment / lagged pass-through: banks close ~55% of the gap
    # to the new long-run target each quarter (Calvo-style repricing friction)
    target_rate = RRP + base_spread
  )

# Partial adjustment recursion for lending rate
pt_df$lending_rate <- NA_real_
pt_df$lending_rate[1] <- pt_df$target_rate[1]
lambda <- 0.55  # speed of adjustment
for (i in 2:nrow(pt_df)) {
  pt_df$lending_rate[i] <- pt_df$lending_rate[i-1] +
    lambda * (pt_df$target_rate[i] - pt_df$lending_rate[i-1]) +
    rnorm(1, 0, 0.05)
}

pt_df <- pt_df %>%
  mutate(
    pass_through_gap = target_rate - lending_rate,   # incomplete PT measure
    d_lending_rate = lending_rate - lag(lending_rate)
  ) %>%
  select(date, t, RRP, d_RRP, lending_rate, d_lending_rate, pass_through_gap)

# -----------------------------------------------------------------------------
# 5. PANEL SKELETON: 2 regions x 4 property types
# -----------------------------------------------------------------------------
regions <- c("NCR", "AONCR")
prop_types <- c("Single-Detached", "Duplex", "Townhouse", "Condominium")

panel_skeleton <- expand_grid(date = quarters_seq, region = regions, property_type = prop_types) %>%
  left_join(time_df, by = "date") %>%
  left_join(macro_df, by = c("date", "t")) %>%
  left_join(pt_df, by = c("date", "t"))

# -----------------------------------------------------------------------------
# 6. STRUCTURAL RATE-SENSITIVITY PARAMETERS BY SEGMENT
#    Priors reflect known institutional facts:
#      - Condos: highest investor/speculative share, most rate-elastic
#      - Single-Detached: owner-occupier dominated, most rate-inelastic
#      - NCR condos specifically: highest financialization, highest beta
#    These become the "true" DGP parameters our Module 2 models should recover
# -----------------------------------------------------------------------------
sensitivity_params <- tribble(
  ~region,  ~property_type,      ~beta_rrp,  ~base_qoq_drift, ~idio_sd,
  "NCR",    "Single-Detached",   -0.55,       0.85,            1.05,
  "NCR",    "Duplex",            -0.70,       0.80,            1.15,
  "NCR",    "Townhouse",         -0.85,       0.90,            1.20,
  "NCR",    "Condominium",       -1.35,       1.00,            1.60,
  "AONCR",  "Single-Detached",   -0.35,       0.75,            0.95,
  "AONCR",  "Duplex",            -0.45,       0.70,            1.00,
  "AONCR",  "Townhouse",         -0.55,       0.78,            1.05,
  "AONCR",  "Condominium",       -0.95,       0.95,            1.35
)

panel <- panel_skeleton %>%
  left_join(sensitivity_params, by = c("region", "property_type")) %>%
  arrange(region, property_type, t) %>%
  group_by(region, property_type) %>%
  mutate(
    # REL (real estate loan) growth: procyclical, dampened by lending rate
    REL_growth = 8 + 0.9 * real_gdp_growth - 1.8 * lending_rate +
      rnorm(n(), 0, 3),
    # Structural QoQ RPPI log-return DGP:
    #   base drift + beta_rrp * lagged rate change (transmission w/ 2-qtr lag)
    #   + macro pass-through (GDP procyclical, inflation weakly positive
    #     via nominal asset repricing) + AR(1) idiosyncratic persistence
    lag1_d_RRP = lag(d_RRP, 1),
    lag2_d_RRP = lag(d_RRP, 2),
    lag4_d_RRP = lag(d_RRP, 4),
    shock = rnorm(n(), 0, idio_sd),
  ) %>%
  mutate(
    dln_RPPI = (base_qoq_drift +
                  beta_rrp * replace_na(lag1_d_RRP, 0) +
                  0.5 * beta_rrp * replace_na(lag2_d_RRP, 0) +
                  0.15 * real_gdp_growth / 4 +
                  0.03 * cpi_national +
                  shock) / 100
  ) %>%
  # Build index levels from log-returns, base = 100 at 2018Q1
  mutate(
    RPPI = 100 * exp(cumsum(replace_na(dln_RPPI, 0)))
  ) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 7. FINALIZE TIDY PANEL + LAGGED TRANSMISSION VARIABLES (t-1, t-2, t-4)
# -----------------------------------------------------------------------------
panel_final <- panel %>%
  group_by(region, property_type) %>%
  arrange(t, .by_group = TRUE) %>%
  mutate(
    RRP_lag1 = lag(RRP, 1), RRP_lag2 = lag(RRP, 2), RRP_lag4 = lag(RRP, 4),
    d_RRP_lag1 = lag(d_RRP, 1), d_RRP_lag2 = lag(d_RRP, 2), d_RRP_lag4 = lag(d_RRP, 4),
    lending_rate_lag1 = lag(lending_rate, 1),
    RPPI_lag1 = lag(RPPI, 1),
    unit_id = paste(region, property_type, sep = "_")
  ) %>%
  ungroup() %>%
  select(
    unit_id, region, property_type, date, yq_lab = t, year, quarter, t,
    RPPI, dln_RPPI, RPPI_lag1,
    RRP, d_RRP, RRP_lag1, RRP_lag2, RRP_lag4,
    d_RRP_lag1, d_RRP_lag2, d_RRP_lag4,
    lending_rate, d_lending_rate, lending_rate_lag1, pass_through_gap,
    REL_growth, cpi_national, real_gdp_growth
  ) %>%
  left_join(time_df %>% select(t, yq_lab_true = yq_lab), by = "t") %>%
  mutate(yq_lab = yq_lab_true) %>%
  select(-yq_lab_true)

# -----------------------------------------------------------------------------
# 8. EXPORT
# -----------------------------------------------------------------------------
dir.create("data", showWarnings = FALSE)
write_csv(panel_final, "data/bsp_rppi_panel.csv")

cat("=== Data generation complete ===\n")
cat("Panel dimensions: N =", n_distinct(panel_final$unit_id),
    " T =", n_distinct(panel_final$t),
    " Total obs =", nrow(panel_final), "\n")
cat("Date range:", as.character(min(panel_final$date)), "to",
    as.character(max(panel_final$date)), "\n")
cat("Saved to data/bsp_rppi_panel.csv\n")
