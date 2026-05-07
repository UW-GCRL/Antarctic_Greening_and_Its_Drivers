"""Stage 2 (regional): Exact Shapley attribution for each region (Continental,
East, West, Peninsula). Reports baseline trend, all-8-detrended trend,
climate-explained fraction, per-variable Shapley, and Mann-Kendall p-values
for baseline vs all-8 for each region. Saves per-region Shapley CSVs."""

from itertools import combinations
from math import factorial

import numpy as np
import pandas as pd
import pymannkendall as mk

from paths import OUTPUT_DIR

VARS = ["temperature_2m", "uv_radiation", "solar_radiation",
        "volumetric_soil_water", "runoff", "X10m_wind_speed",
        "precipitation", "vapor_pressure_deficit"]
N = len(VARS)
ALL = ",".join(sorted(VARS))


def slope_ols(yrs, vals):
    a, _ = np.polyfit(yrs.astype(float), vals, 1)
    return float(a)


def scenario_key(subset):
    return "baseline" if not subset else ",".join(sorted(subset))


def shapley_per_fold(per_fold_csv):
    pf = pd.read_csv(per_fold_csv)
    fold_phis = {var: [] for var in VARS}
    fold_total = []
    fold_baseline = []
    fold_all = []
    for fold in sorted(pf.fold.unique()):
        sub = pf[pf.fold == fold]
        sl = {}
        for sname, g in sub.groupby("scenario"):
            g = g.sort_values("year")
            sl[sname] = slope_ols(g.year.values, g.area.values)
        b = sl["baseline"]

        def v(S):
            return b - sl[scenario_key(S)]

        phi = {var: 0.0 for var in VARS}
        for var in VARS:
            others = [w for w in VARS if w != var]
            for r in range(N):
                wt = factorial(r) * factorial(N - 1 - r) / factorial(N)
                for combo in combinations(others, r):
                    phi[var] += wt * (v(list(combo) + [var]) - v(list(combo)))
        for var in VARS:
            fold_phis[var].append(phi[var])
        fold_total.append(sum(phi.values()))
        fold_baseline.append(b)
        fold_all.append(sl[ALL])
    return fold_phis, fold_total, fold_baseline, fold_all


def mk_significance(per_fold_csv, scenarios=None):
    pf = pd.read_csv(per_fold_csv)
    out = {}
    for sname, g in pf.groupby("scenario"):
        if scenarios is not None and sname not in scenarios:
            continue
        ps = []
        for fold, gg in g.groupby("fold"):
            gg = gg.sort_values("year")
            ps.append(mk.original_test(gg.area.values).p)
        out[sname] = np.mean(ps)
    return out


def main():
    for region in ["Continental", "East", "West", "Peninsula"]:
        pf_path = OUTPUT_DIR / f"{region}_per_fold.csv"
        if not pf_path.exists():
            print(f"\n=== {region}: not yet available ===")
            continue
        print(f"\n========= {region} =========")
        fold_phis, fold_total, fold_baseline, fold_all = shapley_per_fold(pf_path)
        bs = np.mean(fold_baseline)
        a8 = np.mean(fold_all)
        tot = np.mean(fold_total)
        print(f"Baseline slope:           {bs:>8.2f} +/- {np.std(fold_baseline):.2f} km^2/yr")
        print(f"All-8-detrended slope:    {a8:>8.2f} +/- {np.std(fold_all):.2f} km^2/yr")
        print(f"Climate-explained:        {tot:>8.2f} km^2/yr ({100 * tot / bs:.1f}% of baseline trend)")
        print()
        print(f"{'variable':<25} {'phi (km^2/yr)':>16} {'std':>8} {'% of climate-trend':>20}")
        print("-" * 75)
        rows = sorted([(v, np.mean(fold_phis[v]), np.std(fold_phis[v])) for v in VARS],
                      key=lambda r: -r[1])
        for v, m, s in rows:
            pct = 100 * m / tot if tot != 0 else 0
            print(f"{v:<25} {m:>16.2f} {s:>8.2f} {pct:>19.1f}%")

        sigs = mk_significance(pf_path, scenarios=["baseline", ALL])
        print()
        print(f"MK p-value (baseline):    {sigs['baseline']:.4f}  ({'sig' if sigs['baseline'] < 0.05 else 'NOT sig'})")
        print(f"MK p-value (all-8):       {sigs[ALL]:.4f}  ({'sig' if sigs[ALL] < 0.05 else 'NOT sig'})")

        out = pd.DataFrame([{
            "region": region, "variable": v,
            "shapley_mean": np.mean(fold_phis[v]),
            "shapley_std":  np.std(fold_phis[v]),
            "percent_of_climate_trend": 100 * np.mean(fold_phis[v]) / tot if tot != 0 else 0,
        } for v in VARS]).sort_values("shapley_mean", ascending=False)
        out.to_csv(OUTPUT_DIR / f"{region}_shapley.csv", index=False)


if __name__ == "__main__":
    main()
