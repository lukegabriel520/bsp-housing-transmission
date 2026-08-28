#===============================================================================
# 03_model_b_bayesian.R
# MODEL B: Bayesian Hierarchical Model of Property-Specific Rate Sensitivity
#
# SPECIFICATION:
#   dln_RPPI_it ~ Student_t(nu, mu_it, sigma)
#   mu_it = alpha_i + beta_i * d_RRP_lag2_t + gamma * X_it
#   alpha_i ~ Normal(alpha_bar, tau_alpha)          [partial pooling on levels]
#   beta_i  ~ Normal(beta_bar, tau_beta)            [partial pooling on slopes]
#   beta_bar ~ Normal(-0.5, 0.5)                    [weakly informative prior:
#                                                      centered on a modest
#                                                      negative transmission
#                                                      effect, consistent with
#                                                      textbook monetary theory]
#   tau_alpha, tau_beta ~ Student_t+(3, 0, 1)        [half-t hyperpriors]
#   sigma ~ Student_t+(3, 0, 1)
#   nu ~ Gamma(2, 0.1)                                [robust to outlier quarters,
#                                                       e.g. COVID / 2023 tightening]
#
# WHY HIERARCHICAL: With only 8 units x ~30 quarters, unpooled per-unit OLS
#   would be noisy (few effective df per slope). Partial pooling shrinks
#   extreme/noisy unit-level slopes toward the group mean, which is exactly
#   the right prior when we believe all 8 segments share a common monetary
#   transmission mechanism but differ in degree (a classic exchangeability
#   assumption -- property types are functionally similar assets differing
#   mainly in buyer composition / financing structure, not in kind).
#
# OUTPUT: exact 95% posterior credible intervals for each of 8 unit-level
#   beta_i (property-type x region rate sensitivities), plus the population
#   mean effect beta_bar.
#===============================================================================

library(tidyverse)
library(brms)

panel <- read_csv("data/bsp_rppi_panel.csv", show_col_types = FALSE) %>%
  filter(!is.na(dln_RPPI), !is.na(d_RRP_lag2)) %>%
  mutate(
    unit_id = factor(unit_id),
    property_type = factor(property_type,
                            levels = c("Single-Detached", "Duplex", "Townhouse", "Condominium")),
    region = factor(region),
    # standardize macro controls for sane prior scales / faster sampling
    cpi_z = as.numeric(scale(cpi_national)),
    gdp_z = as.numeric(scale(real_gdp_growth))
  )

# -----------------------------------------------------------------------------
# PRIORS
#   Student-t(3, 0, 1) half-priors on all variance/hyperparameters are a
#   standard robust "weakly informative" choice (Gelman 2006) that regularizes
#   without dominating the likelihood given our N=8, T~30 data.
# -----------------------------------------------------------------------------
priors <- c(
  prior(normal(-0.5, 0.5), class = "b", coef = "d_RRP_lag2"),
  prior(normal(0, 1),      class = "b", coef = "cpi_z"),
  prior(normal(0, 1),      class = "b", coef = "gdp_z"),
  prior(student_t(3, 0, 1), class = "sd"),                # random-effect SDs
  prior(student_t(3, 0, 1), class = "sigma"),              # residual SD
  prior(gamma(2, 0.1),      class = "nu")                  # Student-t df (robustness)
)

# -----------------------------------------------------------------------------
# MODEL: varying intercept AND varying slope on d_RRP_lag2 by unit_id
#   (random-effects hierarchical / "multilevel" model)
# -----------------------------------------------------------------------------
mod_b <- brm(
  bf(dln_RPPI ~ d_RRP_lag2 + cpi_z + gdp_z + (1 + d_RRP_lag2 | unit_id)),
  data = panel,
  family = student(),           # Student-t likelihood for robustness to outlier quarters
  prior = priors,
  chains = 4,
  iter = 4000,
  warmup = 1500,
  cores = 4,
  seed = 20260828,
  control = list(adapt_delta = 0.97, max_treedepth = 12),
  backend = "cmdstanr"
)

cat("\n============================================================\n")
cat("MODEL B: Bayesian Hierarchical Model Summary\n")
cat("============================================================\n")
print(summary(mod_b))

# -----------------------------------------------------------------------------
# CONVERGENCE DIAGNOSTICS
# -----------------------------------------------------------------------------
rhat_max <- max(brms::rhat(mod_b), na.rm = TRUE)
ess_min  <- min(brms::neff_ratio(mod_b) * (4000 - 1500) * 4, na.rm = TRUE)
cat("\nMax Rhat:", round(rhat_max, 4), " (should be < 1.01)\n")
cat("Min effective sample size:", round(ess_min, 0), " (should be > 400)\n")

# -----------------------------------------------------------------------------
# EXTRACT UNIT-LEVEL POSTERIOR CREDIBLE INTERVALS (exact 95% CrI via quantiles
# of the posterior draws -- not a normal-approximation interval)
# -----------------------------------------------------------------------------
unit_lookup <- panel %>% distinct(unit_id, region, property_type)

ranef_draws <- ranef(mod_b, summary = FALSE)$unit_id[, , "d_RRP_lag2"]
fixef_draws <- as_draws_df(mod_b)$b_d_RRP_lag2

unit_beta_draws <- sweep(ranef_draws, 1, fixef_draws, "+")  # total slope per unit per draw

posterior_summary <- tibble(
  unit_id = colnames(unit_beta_draws)
) %>%
  left_join(unit_lookup, by = "unit_id") %>%
  mutate(
    posterior_mean = apply(unit_beta_draws, 2, mean),
    posterior_median = apply(unit_beta_draws, 2, median),
    cri_2.5  = apply(unit_beta_draws, 2, quantile, probs = 0.025),
    cri_97.5 = apply(unit_beta_draws, 2, quantile, probs = 0.975),
    prob_negative = apply(unit_beta_draws, 2, function(x) mean(x < 0))
  ) %>%
  arrange(posterior_mean)

cat("\n============================================================\n")
cat("Unit-level posterior rate-sensitivity (beta_i), exact 95% credible intervals\n")
cat("============================================================\n")
print(posterior_summary, n = Inf)

write_csv(posterior_summary, "data/bayesian_posterior_summary.csv")
saveRDS(mod_b, "data/model_b_object.rds")

# -----------------------------------------------------------------------------
# POPULATION-LEVEL (POOLED) EFFECT
# -----------------------------------------------------------------------------
pop_effect <- as_draws_df(mod_b)$b_d_RRP_lag2
cat("\nPopulation mean effect (beta_bar):\n")
cat("  Posterior mean:", round(mean(pop_effect), 4), "\n")
cat("  95% CrI: [", round(quantile(pop_effect, 0.025), 4), ",",
    round(quantile(pop_effect, 0.975), 4), "]\n")
cat("  P(beta_bar < 0):", round(mean(pop_effect < 0), 4), "\n")

cat("\nModel B complete. Posterior summaries saved to data/bayesian_posterior_summary.csv\n")
