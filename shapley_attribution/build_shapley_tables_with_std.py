"""Build region x variable Shapley tables with mean +/- std across 5 folds.
Saves two CSVs: absolute (km^2/yr) and percentage (% of climate-explained)."""

from itertools import combinations
from math import factorial

import numpy as np
import pandas as pd

from paths import OUTPUT_DIR

VARS_DISPLAY = ["volumetric_soil_water", "precipitation", "uv_radiation", "temperature_2m",
                "X10m_wind_speed", "solar_radiation", "runoff", "vapor_pressure_deficit"]
VARS = VARS_DISPLAY[:]
N = len(VARS)
SHORT = {"volumetric_soil_water": "VSW", "precipitation": "Precip", "uv_radiation": "UV",
         "temperature_2m": "T_2m", "X10m_wind_speed": "Wind", "solar_radiation": "Solar",
         "runoff": "Runoff", "vapor_pressure_deficit": "VPD"}
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
    return fold_phis, fold_total


def main():
    regions = ["Continental", "East", "West", "Peninsula"]
    all_results = {}
    for r in regions:
        fold_phis, fold_total = shapley_per_fold(OUTPUT_DIR / f"{r}_per_fold.csv")
        all_results[r] = {"phis": fold_phis, "total": fold_total}

    cols = [SHORT[v] for v in VARS_DISPLAY]

    # Absolute
    print("=== Shapley attribution (km^2/yr, mean +/- std across 5 folds) ===")
    header = f"{'region':<12} | " + " | ".join(f"{c:>15}" for c in cols) + f" | {'TOTAL':>15}"
    print(header)
    print("-" * len(header))
    for r in regions:
        means = [np.mean(all_results[r]["phis"][v]) for v in VARS_DISPLAY]
        stds  = [np.std(all_results[r]["phis"][v])  for v in VARS_DISPLAY]
        total_m = np.mean(all_results[r]["total"])
        total_s = np.std(all_results[r]["total"])
        cells = [f"{m:+6.2f} +/- {s:5.2f}" for m, s in zip(means, stds)]
        print(f"{r:<12} | " + " | ".join(f"{c:>15}" for c in cells)
              + f" | {f'{total_m:+6.2f} +/- {total_s:5.2f}':>15}")

    # Percentage
    print()
    print("=== Shapley attribution (% of region's climate-explained trend) ===")
    print(header)
    print("-" * len(header))
    for r in regions:
        total_per_fold = np.array(all_results[r]["total"])
        rows_pct = []
        for v in VARS_DISPLAY:
            per_fold_pct = 100 * np.array(all_results[r]["phis"][v]) / total_per_fold
            rows_pct.append((per_fold_pct.mean(), per_fold_pct.std()))
        cells = [f"{m:+6.1f}% +/-{s:4.1f}" for m, s in rows_pct]
        print(f"{r:<12} | " + " | ".join(f"{c:>15}" for c in cells) + f" | {'100.0%':>15}")

    # Save CSVs
    abs_rows, pct_rows = [], []
    for r in regions:
        means = {SHORT[v]: np.mean(all_results[r]["phis"][v]) for v in VARS_DISPLAY}
        stds  = {SHORT[v]: np.std(all_results[r]["phis"][v]) for v in VARS_DISPLAY}
        total_m = np.mean(all_results[r]["total"])
        total_s = np.std(all_results[r]["total"])
        abs_row = {"region": r}
        for c in cols:
            abs_row[f"{c}_mean"] = means[c]
            abs_row[f"{c}_std"] = stds[c]
        abs_row["TOTAL_mean"] = total_m
        abs_row["TOTAL_std"] = total_s
        abs_rows.append(abs_row)

        total_per_fold = np.array(all_results[r]["total"])
        pct_row = {"region": r}
        for v in VARS_DISPLAY:
            per_fold_pct = 100 * np.array(all_results[r]["phis"][v]) / total_per_fold
            pct_row[f"{SHORT[v]}_mean"] = per_fold_pct.mean()
            pct_row[f"{SHORT[v]}_std"] = per_fold_pct.std()
        pct_rows.append(pct_row)

    pd.DataFrame(abs_rows).to_csv(OUTPUT_DIR / "shapley_table_abs_with_std.csv", index=False)
    pd.DataFrame(pct_rows).to_csv(OUTPUT_DIR / "shapley_table_pct_with_std.csv", index=False)
    print()
    print(f"Saved: shapley_table_abs_with_std.csv  and  shapley_table_pct_with_std.csv")


if __name__ == "__main__":
    main()
