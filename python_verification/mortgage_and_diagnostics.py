"""
mortgage_and_diagnostics.py
Computes the mortgage stress matrix (exact amortization math, no simulation
needed) and runs verification diagnostics on the synthetic panel.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.diagnostic import het_breuschpagan
from statsmodels.tsa.stattools import adfuller
import statsmodels.formula.api as smf
import statsmodels.api as sm

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
DATA.mkdir(parents=True, exist_ok=True)

panel = pd.read_csv(DATA / "bsp_rppi_panel.csv")

# =============================================================================
# MORTGAGE STRESS MATRIX (exact PMT formula, not simulated)
# =============================================================================
principal = 5_000_000
years = 20
n_months = years * 12

latest_t = panel["t"].max()
baseline_mortgage_rate = panel.loc[panel["t"] == latest_t, "lending_rate"].mean()

def monthly_amortization(principal, annual_rate_pct, n_months):
    r = (annual_rate_pct / 100) / 12
    if r == 0:
        return principal / n_months
    return principal * r * (1 + r) ** n_months / ((1 + r) ** n_months - 1)

scenarios = [("Baseline", 0), ("+50 bps", 50), ("+100 bps", 100), ("+200 bps", 200)]
rows = []
baseline_amort = None
for name, bps in scenarios:
    rate = baseline_mortgage_rate + bps / 100
    amort = monthly_amortization(principal, rate, n_months)
    if name == "Baseline":
        baseline_amort = amort
    total_interest = amort * n_months - principal
    rows.append({
        "scenario": name, "mortgage_rate": rate, "monthly_amortization": amort,
        "total_interest_paid": total_interest
    })

mort_df = pd.DataFrame(rows)
mort_df["delta_vs_baseline"] = mort_df["monthly_amortization"] - baseline_amort
mort_df["pct_increase_vs_baseline"] = 100 * mort_df["delta_vs_baseline"] / baseline_amort

print("=== MORTGAGE STRESS MATRIX ===")
print(f"Baseline mortgage rate anchor: {baseline_mortgage_rate:.2f}%")
print(mort_df.round(2).to_string(index=False))
mort_df.to_csv(DATA / "mortgage_stress_matrix_verified.csv", index=False)

# =============================================================================
# DIAGNOSTIC TESTS
# =============================================================================
print("\n=== 1. AUGMENTED DICKEY-FULLER (panel-by-panel, proxy for LLC/IPS logic) ===")
# Run ADF per unit on ln(RPPI) level and on dln_RPPI, report % rejecting unit root
adf_level_pvals = []
adf_diff_pvals = []
for uid, grp in panel.groupby("unit_id"):
    grp = grp.sort_values("t")
    ln_rppi = np.log(grp["RPPI"].values)
    try:
        adf_level_pvals.append(adfuller(ln_rppi, autolag="AIC")[1])
    except Exception:
        pass
    dln = grp["dln_RPPI"].dropna().values
    try:
        adf_diff_pvals.append(adfuller(dln, autolag="AIC")[1])
    except Exception:
        pass

print(f"ln(RPPI) levels: {sum(p < 0.05 for p in adf_level_pvals)}/{len(adf_level_pvals)} units reject unit-root null at 5% (expect: FEW -- confirms I(1))")
print(f"  mean p-value: {np.mean(adf_level_pvals):.4f}")
print(f"Delta ln(RPPI): {sum(p < 0.05 for p in adf_diff_pvals)}/{len(adf_diff_pvals)} units reject unit-root null at 5% (expect: MOST/ALL -- confirms I(0))")
print(f"  mean p-value: {np.mean(adf_diff_pvals):.6f}")

print("\n=== 2. SERIAL CORRELATION: Breusch-Godfrey-style check via residual autocorrelation ===")
reg_data = panel.dropna(subset=["dln_RPPI", "d_RRP_lag2"]).copy()
m = smf.ols("dln_RPPI ~ d_RRP_lag2 + cpi_national + real_gdp_growth + C(unit_id)", data=reg_data).fit()
resid = m.resid
# unit-by-unit lag-1 autocorrelation of residuals
autocorrs = []
for uid, idx in reg_data.groupby("unit_id").groups.items():
    r = resid.loc[idx].values
    if len(r) > 3:
        ac1 = np.corrcoef(r[:-1], r[1:])[0, 1]
        autocorrs.append(ac1)
print(f"Mean within-unit lag-1 residual autocorrelation: {np.mean(autocorrs):.4f}")
print("(Nonzero autocorrelation supports using cluster-robust / wild-bootstrap SEs)")

print("\n=== 3. HETEROSKEDASTICITY: Breusch-Pagan ===")
X = sm.add_constant(reg_data[["d_RRP_lag2", "cpi_national", "real_gdp_growth"]])
bp_stat, bp_pval, f_stat, f_pval = het_breuschpagan(m.resid, X)
print(f"BP LM statistic: {bp_stat:.3f}, p-value: {bp_pval:.4f}")
print(f"BP F statistic: {f_stat:.3f}, p-value: {f_pval:.4f}")

print("\n=== 4. PASS-THROUGH COMPLETENESS: RRP -> Lending Rate ===")
pt_data = panel.drop_duplicates(subset=["t"])[["t", "RRP", "lending_rate", "d_RRP", "d_lending_rate"]].dropna()
pt_model = smf.ols("d_lending_rate ~ d_RRP", data=pt_data).fit()
print(pt_model.summary().tables[1])
pt_coef = pt_model.params["d_RRP"]
pt_se = pt_model.bse["d_RRP"]
t_complete = (pt_coef - 1) / pt_se
from scipy import stats as sstats
p_complete = 2 * (1 - sstats.t.cdf(abs(t_complete), df=pt_model.df_resid))
print(f"\nContemporaneous pass-through coefficient: {pt_coef:.3f}")
print(f"Test of H0 (complete pass-through, coef=1): t={t_complete:.3f}, p={p_complete:.4f}")
print("=> " + ("REJECT: pass-through incomplete same-quarter" if p_complete < 0.05 else "FAIL TO REJECT: consistent with complete pass-through"))

# 4-quarter cumulative pass-through
pt_lr = panel.drop_duplicates(subset=["t"])[["t", "RRP", "lending_rate"]].sort_values("t").copy()
pt_lr["RRP_chg_4q"] = pt_lr["RRP"].diff(4)
pt_lr["lend_chg_4q"] = pt_lr["lending_rate"].diff(4)
pt_lr = pt_lr.dropna()
pt_lr_model = smf.ols("lend_chg_4q ~ RRP_chg_4q", data=pt_lr).fit()
print(f"\nLong-run (4Q) pass-through coefficient: {pt_lr_model.params['RRP_chg_4q']:.3f}")
print(f"  (partial-adjustment lambda=0.55 by construction implies theoretical "
      f"4Q cumulative pass-through of ~{1-(1-0.55)**4:.3f} of the long-run spread target)")

diag_summary = pd.DataFrame({
    "test": ["ADF level (mean p)", "ADF diff (mean p)", "Mean residual autocorr",
             "Breusch-Pagan LM", "Contemporaneous pass-through", "Long-run 4Q pass-through"],
    "statistic": [np.mean(adf_level_pvals), np.mean(adf_diff_pvals), np.mean(autocorrs),
                  bp_stat, pt_coef, pt_lr_model.params["RRP_chg_4q"]],
    "p_value": [np.nan, np.nan, np.nan, bp_pval, pt_model.pvalues["d_RRP"], pt_lr_model.pvalues["RRP_chg_4q"]]
})
diag_summary.to_csv(DATA / "diagnostic_summary_verified.csv", index=False)
print("\nSaved diagnostic_summary_verified.csv")
