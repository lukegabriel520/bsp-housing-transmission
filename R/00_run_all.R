#===============================================================================
# 00_run_all.R
# MASTER PIPELINE RUNNER
# BSP Monetary Policy -> Philippine Residential Real Estate Transmission Study
#
# Run this single script to execute the full pipeline end-to-end:
#   1. Generate calibrated synthetic panel data
#   2. Estimate Model A (frequentist FE panel + wild cluster bootstrap + LP-IRF)
#   3. Estimate Model B (Bayesian hierarchical model)
#   4. Run diagnostics (unit root, serial correlation, heteroskedasticity,
#      pass-through completeness)
#   5. Generate IRF plots and mortgage stress matrix
#
# REQUIRED PACKAGES:
#   install.packages(c("tidyverse","fixest","fwildclusterboot","modelsummary",
#                       "brms","plm","lmtest","patchwork","scales","lubridate"))
#   # brms requires a working Stan backend; recommended:
#   install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
#   cmdstanr::install_cmdstan()
#===============================================================================

script_path <- tryCatch(
  rstudioapi::getActiveDocumentContext()$path,
  error = function(e) NA_character_
)
if (!is.na(script_path) && nzchar(script_path)) {
  setwd(normalizePath(file.path(dirname(script_path), "..")))
} else if (file.exists("R/00_run_all.R")) {
  # already at project root when run via Rscript from root
} else if (file.exists("00_run_all.R")) {
  setwd(normalizePath(file.path(getwd(), "..")))
}
# setwd("/path/to/bsp-housing-transmission")  # <-- uncomment and edit if needed

dir.create("data", showWarnings = FALSE)
dir.create("outputs", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)
dir.create("report", showWarnings = FALSE)

cat("\n#####################################################################\n")
cat("# STEP 1/5: DATA GENERATION\n")
cat("#####################################################################\n")
source("R/01_generate_data.R")

cat("\n#####################################################################\n")
cat("# STEP 2/5: MODEL A (Frequentist Panel + Wild Cluster Bootstrap + LP-IRF)\n")
cat("#####################################################################\n")
source("R/02_model_a_frequentist.R")

cat("\n#####################################################################\n")
cat("# STEP 3/5: MODEL B (Bayesian Hierarchical Model)\n")
cat("#####################################################################\n")
source("R/03_model_b_bayesian.R")

cat("\n#####################################################################\n")
cat("# STEP 4/5: DIAGNOSTICS\n")
cat("#####################################################################\n")
source("R/04_diagnostics.R")

cat("\n#####################################################################\n")
cat("# STEP 5/5: IRF PLOTS & MORTGAGE STRESS MATRIX\n")
cat("#####################################################################\n")
source("R/05_irf_plot_and_mortgage_stress.R")

cat("\n#####################################################################\n")
cat("# PIPELINE COMPLETE\n")
cat("#####################################################################\n")
cat("Outputs:\n")
cat("  data/bsp_rppi_panel.csv               - full tidy panel dataset\n")
cat("  data/wild_cluster_bootstrap_results.csv\n")
cat("  data/local_projections_irf.csv\n")
cat("  data/regression_results.csv\n")
cat("  data/bayesian_posterior_summary.csv\n")
cat("  data/diagnostic_tests.csv\n")
cat("  data/mortgage_stress_matrix.csv\n")
cat("  outputs/table1_regression_results.docx\n")
cat("  outputs/mortgage_stress_matrix.md\n")
cat("  plots/irf_faceted.png\n")
cat("  plots/irf_overlay.png\n")
cat("  plots/mortgage_stress_matrix.png\n")
