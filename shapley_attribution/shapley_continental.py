"""Stage 2 (continental): Exact Shapley value attribution for the 8 climate
variables on the continental greening trend.

Value function:
    v(S) = baseline_slope - slope(scenario where vars in S are detrended)
That is, "amount of trend explained by removing the trends of variables in S".
    v({})       = 0
    v({all 8})  = baseline_slope - slope(all_8) = total trend explained
    Σ_v φ_v     = v({all 8})  by Shapley efficiency.

Uncertainty: compute Shapley per-fold, then report mean ± std across folds.
"""

from itertools import combinations
from math import factorial

import numpy as np
import pandas as pd

from paths import OUTPUT_DIR

VARS = ["temperature_2m", "uv_radiation", "solar_radiation",
        "volumetric_soil_water", "runoff", "X10m_wind_speed",
        "precipitation", "vapor_pressure_deficit"]
N = len(VARS)


def slope_ols(yrs, vals):
    a, _ = np.polyfit(yrs.astype(float), vals, 1)
    return float(a)


def scenario_key(subset):
    if not subset:
        return "baseline"
    return ",".join(sorted(subset))


def compute_value_lookup_for_fold(pf, fold):
    """Returns dict scenario_name -> OLS slope of yearly area for that scenario in this fold."""
    sub = pf[pf.fold == fold]
    out = {}
    for sname, g in sub.groupby("scenario"):
        g = g.sort_values("year")
        out[sname] = slope_ols(g.year.values, g.area.values)
    return out


def shapley(slopes_lookup):
    """Given a dict scenario_name -> slope, compute exact Shapley per variable."""
    baseline_slope = slopes_lookup["baseline"]

    def v(S):
        return baseline_slope - slopes_lookup[scenario_key(S)]

    phi = {var: 0.0 for var in VARS}
    for var in VARS:
        others = [w for w in VARS if w != var]
        for r in range(N):  # subset size 0..N-1
            wt = factorial(r) * factorial(N - 1 - r) / factorial(N)
            for combo in combinations(others, r):
                S = list(combo)
                S_with = list(combo) + [var]
                phi[var] += wt * (v(S_with) - v(S))
    return phi, baseline_slope


def main():
    pf = pd.read_csv(OUTPUT_DIR / "Continental_per_fold.csv")

    fold_phis = {var: [] for var in VARS}
    fold_total = []
    fold_baselines = []
    fold_all_detrended = []
    all_key = ",".join(sorted(VARS))
    for fold in sorted(pf.fold.unique()):
        sl = compute_value_lookup_for_fold(pf, fold)
        phi, b = shapley(sl)
        for var in VARS:
            fold_phis[var].append(phi[var])
        fold_total.append(sum(phi.values()))
        fold_baselines.append(b)
        fold_all_detrended.append(sl[all_key])

    print("Continental Shapley attribution of greening trend (km^2/yr)")
    print("Value function: trend reduction when variable's trend is removed.")
    print("Mean +/- std across 5 folds.")
    print()
    print(f"{'variable':<25} {'phi (km^2/yr)':>16} {'std':>10} {'% of total':>12}")
    print("-" * 70)
    mean_total = np.mean(fold_total)
    rows = []
    for var in VARS:
        arr = np.array(fold_phis[var])
        rows.append((var, arr.mean(), arr.std(), 100 * arr.mean() / mean_total))
    rows.sort(key=lambda r: -r[1])
    for var, m, s, pct in rows:
        print(f"{var:<25} {m:>16.3f} {s:>10.3f} {pct:>11.1f}%")
    print("-" * 70)
    print(f"{'TOTAL (sum of phi)':<25} {mean_total:>16.3f} {np.std(fold_total):>10.3f} {100.0:>11.1f}%")
    print()
    print(f"Baseline slope (mean across folds):           {np.mean(fold_baselines):.2f} km^2/yr")
    print(f"All-8-detrended slope (mean across folds):    {np.mean(fold_all_detrended):.2f} km^2/yr")
    print(f"Total trend explained (baseline - all8):      "
          f"{np.mean(fold_baselines) - np.mean(fold_all_detrended):.2f} km^2/yr")
    print("Shapley sum should match the line above (efficiency property).")

    out = pd.DataFrame([{
        "variable": var,
        "shapley_mean_km2_per_yr": np.mean(fold_phis[var]),
        "shapley_std_km2_per_yr":  np.std(fold_phis[var]),
        "percent_of_total":        100 * np.mean(fold_phis[var]) / mean_total,
    } for var in VARS]).sort_values("shapley_mean_km2_per_yr", ascending=False)
    out_path = OUTPUT_DIR / "Continental_shapley.csv"
    out.to_csv(out_path, index=False)
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
