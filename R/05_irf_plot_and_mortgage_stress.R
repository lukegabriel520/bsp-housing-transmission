#===============================================================================
# 05_irf_plot_and_mortgage_stress.R
# MODULE 4 (viz) + MODULE 5.2: IRF Plot & Mortgage Stress Matrix
#===============================================================================

library(tidyverse)
library(patchwork)
library(scales)

# -----------------------------------------------------------------------------
# A. IMPULSE RESPONSE FUNCTION PLOT
#    8-quarter price trajectory of each property type following a +100bp hike
# -----------------------------------------------------------------------------
lp_results <- read_csv("data/local_projections_irf.csv", show_col_types = FALSE)

irf_colors <- c(
  "Single-Detached" = "#2C5F2D",
  "Duplex"          = "#5B8C5A",
  "Townhouse"        = "#E8A33D",
  "Condominium"      = "#C1440E"
)

p_irf <- ggplot(lp_results, aes(x = h, y = irf_pct, color = property_type, fill = property_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_ribbon(aes(ymin = irf_ci_low, ymax = irf_ci_high), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.8) +
  scale_color_manual(values = irf_colors) +
  scale_fill_manual(values = irf_colors) +
  scale_x_continuous(breaks = 0:8) +
  labs(
    title = "Impulse Response: Residential Property Prices to a +100bp BSP Rate Hike",
    subtitle = "Cumulative % change in RPPI by horizon (quarters), local projections, 95% CI",
    x = "Quarters after shock",
    y = "Cumulative RPPI response (%)",
    color = "Property Type", fill = "Property Type",
    caption = "Model: panel local projections, unit (region x type) fixed effects,\ncluster-robust SEs (unit_id, N=8). Shaded bands = 95% CI."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey30", size = 10),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  facet_wrap(~ property_type, nrow = 2)

ggsave("plots/irf_faceted.png", p_irf, width = 10, height = 7, dpi = 300, bg = "white")

# Combined single-panel overlay version too
p_irf_overlay <- ggplot(lp_results, aes(x = h, y = irf_pct, color = property_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  scale_color_manual(values = irf_colors) +
  scale_x_continuous(breaks = 0:8) +
  labs(
    title = "8-Quarter Price Trajectory Following a +100bp BSP Policy Rate Hike",
    subtitle = "By property type, cumulative % RPPI response (point estimates)",
    x = "Quarters after shock", y = "Cumulative RPPI response (%)",
    color = "Property Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("plots/irf_overlay.png", p_irf_overlay, width = 9, height = 6, dpi = 300, bg = "white")

cat("IRF plots saved to plots/irf_faceted.png and plots/irf_overlay.png\n")

# -----------------------------------------------------------------------------
# B. MORTGAGE STRESS MATRIX
#    Standard PHP 5,000,000, 20-year mortgage under +50/+100/+200bp scenarios
# -----------------------------------------------------------------------------
principal <- 5000000
years <- 20
n_months <- years * 12

# Baseline mortgage rate assumption: current typical PH bank housing loan
# fixed rate ~ lending_rate + small retail spread. Use latest lending_rate
# from panel as the anchor.
panel_raw <- read_csv("data/bsp_rppi_panel.csv", show_col_types = FALSE)
baseline_mortgage_rate <- panel_raw %>%
  filter(t == max(t)) %>%
  summarise(r = mean(lending_rate, na.rm = TRUE)) %>%
  pull(r)

cat("\nBaseline mortgage rate anchor (from latest panel lending_rate):",
    round(baseline_mortgage_rate, 2), "%\n")

monthly_amortization <- function(principal, annual_rate_pct, n_months) {
  r <- (annual_rate_pct / 100) / 12
  if (r == 0) return(principal / n_months)
  principal * r * (1 + r)^n_months / ((1 + r)^n_months - 1)
}

scenarios <- tribble(
  ~scenario,        ~bps_change,
  "Baseline",        0,
  "+50 bps",         50,
  "+100 bps",        100,
  "+200 bps",        200
)

mortgage_stress <- scenarios %>%
  mutate(
    mortgage_rate = baseline_mortgage_rate + bps_change / 100,
    monthly_amortization = map_dbl(mortgage_rate, ~monthly_amortization(principal, .x, n_months)),
    total_interest_paid = monthly_amortization * n_months - principal
  ) %>%
  mutate(
    delta_vs_baseline = monthly_amortization - monthly_amortization[scenario == "Baseline"],
    pct_increase_vs_baseline = 100 * delta_vs_baseline / monthly_amortization[scenario == "Baseline"]
  )

print(mortgage_stress)
write_csv(mortgage_stress, "data/mortgage_stress_matrix.csv")

# Format as Markdown table for the civilian report
mortgage_md <- mortgage_stress %>%
  transmute(
    Scenario = scenario,
    `Mortgage Rate` = paste0(round(mortgage_rate, 2), "%"),
    `Monthly Amortization` = paste0("PHP ", format(round(monthly_amortization), big.mark = ",")),
    `Increase vs. Baseline` = ifelse(scenario == "Baseline", "—",
                                      paste0("+PHP ", format(round(delta_vs_baseline), big.mark = ","),
                                             " (", sprintf("%.1f", pct_increase_vs_baseline), "%)")),
    `Total Interest (20yr)` = paste0("PHP ", format(round(total_interest_paid), big.mark = ","))
  )

writeLines(
  c(
    "| Scenario | Mortgage Rate | Monthly Amortization | Increase vs. Baseline | Total Interest (20yr) |",
    "|---|---|---|---|---|",
    apply(mortgage_md, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  ),
  "outputs/mortgage_stress_matrix.md"
)

cat("\nMortgage stress matrix saved to data/mortgage_stress_matrix.csv\n")
cat("and outputs/mortgage_stress_matrix.md\n")

# -----------------------------------------------------------------------------
# C. MORTGAGE STRESS VISUALIZATION
# -----------------------------------------------------------------------------
p_mortgage <- ggplot(mortgage_stress, aes(x = fct_inorder(scenario), y = monthly_amortization)) +
  geom_col(fill = "#2C5F2D", width = 0.6) +
  geom_text(aes(label = paste0("PHP ", format(round(monthly_amortization), big.mark = ","))),
            vjust = -0.5, size = 3.8, fontface = "bold") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Mortgage Stress Matrix: PHP 5,000,000 / 20-Year Loan",
    subtitle = paste0("Monthly amortization under BSP rate-hike scenarios (baseline mortgage rate: ",
                       round(baseline_mortgage_rate, 2), "%)"),
    x = NULL, y = "Monthly Amortization (PHP)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

ggsave("plots/mortgage_stress_matrix.png", p_mortgage, width = 8, height = 5.5, dpi = 300, bg = "white")

cat("Mortgage stress plot saved to plots/mortgage_stress_matrix.png\n")
