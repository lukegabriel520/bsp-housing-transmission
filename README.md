# BSP Housing Transmission Study

Hi there. This repo is my attempt to answer a simple question with real data tools: when the Bangko Sentral ng Pilipinas (BSP) moves interest rates, how long does it take for home prices to feel it, and which types of homes move the most?

I built this because rate hikes kept making headlines, and I wanted to see the link between policy, mortgage payments, and property prices in plain numbers. Condos and single-detached homes do not react the same way. I wanted to show that. (P.S I intend to move out soon enough so renting prices in high-rated cities would extremely empty my alr empty pockets)

## Why this exists

BSP publishes useful summary stats on residential property prices (RPPI), but the loan-level detail behind those numbers is not public at the quarterly region-by-property-type level I needed. So I built a **calibrated synthetic panel**: fake microdata shaped to match real BSP and PSA turning points from 2018 through early 2026. The goal is reproducibility, not a claim that this is official BSP microdata.

The study covers:

- How BSP policy rate (RRP) changes pass through to lending rates and home price growth
- Differences across NCR vs rest of Philippines, and across four property types
- A mortgage stress table: what a rate hike does to a typical monthly payment
- Charts showing price paths after a +100 basis point rate shock (an **impulse response**, or IRF for short)

## What is in each folder

```
bsp-housing-transmission/
├── R/                      # Main analysis pipeline (run this first)
├── python_verification/    # Cross-checks key R numbers in Python
├── data/                   # Panel dataset and result tables (CSV)
├── plots/                  # IRF and mortgage stress charts (PNG)
├── report/                 # Written study (Word doc)
└── outputs/                # Auto-generated tables when you re-run R (not in git)
```

| Folder | What you will find |
|--------|-------------------|
| `R/` | Six scripts, numbered in order. Start with `00_run_all.R` to run everything. |
| `python_verification/` | Optional sanity checks. Writes `*_verified.csv` files into `data/`. |
| `data/` | The panel (`bsp_rppi_panel.csv`) plus regression, bootstrap, IRF, diagnostic, and mortgage outputs. |
| `plots/` | `irf_faceted.png` (by property type), `irf_overlay.png` (all types together), `mortgage_stress_matrix.png`. |
| `report/` | `BSP_Housing_Transmission_Study.docx`, the full writeup. |

## How to run

### R pipeline (full study)

From the project root:

```bash
Rscript R/00_run_all.R
```

Or open `R/00_run_all.R` in RStudio and run it. The script sets the working directory to the project root automatically.

**Required R packages:**

```r
install.packages(c(
  "tidyverse", "fixest", "fwildclusterboot", "modelsummary",
  "brms", "plm", "lmtest", "patchwork", "scales", "lubridate"
))
```

Model B (Bayesian) needs a Stan backend. The script header in `R/00_run_all.R` has `cmdstanr` install notes.

### Python verification (optional)

From the project root:

```bash
python python_verification/verify_pipeline.py
python python_verification/mortgage_and_diagnostics.py
```

These read and write under `data/`. Verified outputs use the `_verified` suffix and are gitignored.

**Python packages:** `numpy`, `pandas`, `statsmodels`, `scipy`

## Key outputs

| File | What it is |
|------|-----------|
| `data/bsp_rppi_panel.csv` | Quarterly panel: 8 region x property-type units, 34 quarters |
| `data/regression_results.csv` | Model A1 interaction coefficients (policy rate x property type) |
| `data/wild_cluster_bootstrap_results.csv` | Small-sample robust p-values (wild cluster bootstrap) |
| `data/local_projections_irf.csv` | 8-quarter price response after a +100bp hike |
| `data/diagnostic_tests.csv` | Unit root, serial correlation, heteroskedasticity, pass-through checks |
| `data/mortgage_stress_matrix.csv` | Monthly payment impact for +0, +50, +100, +200bp scenarios |

## Important note on the data

All panel observations are **synthetic but calibrated** to published BSP RRP paths and PSA/BSP RPPI trends. Use this repo to reproduce the analysis logic and explore transmission patterns. Do not treat the CSV as confidential BSP loan records.

## License

MIT License. See [LICENSE](LICENSE).
