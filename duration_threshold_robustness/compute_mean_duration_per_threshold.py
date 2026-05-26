"""
Robustness check: mean greenness duration trend under different
duration-day thresholds.

For each year (2002-2023), in two duration datasets:
  1. First_last/Land  -> first-to-last NDVI>0.2 day method
  2. all_whole_map/Land -> daily-count method
keep pixels whose duration exceeds each of T = 1..14 days, take the
mean duration across the kept pixels, and run Mann-Kendall +
Sen's slope + linear regression on the 22-year time series.

(All 500 m MODIS pixels are equal-area, so a simple pixel-mean is
equivalent to the area-weighted mean.)

Outputs:
  duration_threshold_mean_results.csv  -- per (dataset, threshold, year)
  duration_threshold_mean_summary.csv  -- per (dataset, threshold)
  mean_duration_by_year_<dataset>.csv  -- pivoted (year x threshold)
"""

from pathlib import Path
import numpy as np
import pandas as pd
import rasterio
import pymannkendall as mk
from scipy import stats

ROOT = Path(r"G:/Hangkai/Anttarctic Vegetation Dynamic/Version_2_data/Antarctica_Vegetation_Duration")

DATASETS = {
    "first_last": {
        "dir":  ROOT / "First_last" / "Land",
        "pattern": "Duration_{year}.tif",
    },
    "daily_count": {
        "dir":  ROOT / "all_whole_map" / "Land",
        "pattern": "Land_{year}.tif",
    },
}

YEARS = list(range(2002, 2024))
THRESHOLDS = list(range(1, 15))   # pixel kept if duration > T days

def mean_duration_for_year(tif_path: Path, t_days: int):
    with rasterio.open(tif_path) as ds:
        arr = ds.read(1)
    mask = arr > t_days
    n = int(mask.sum())
    if n == 0:
        return np.nan, 0
    return float(arr[mask].mean()), n

rows = []
for ds_name, cfg in DATASETS.items():
    for year in YEARS:
        tif = cfg["dir"] / cfg["pattern"].format(year=year)
        if not tif.exists():
            print(f"missing: {tif}")
            continue
        for t in THRESHOLDS:
            mean_dur, npx = mean_duration_for_year(tif, t)
            rows.append({
                "dataset": ds_name,
                "threshold_days": t,
                "year": year,
                "n_pixels": npx,
                "mean_duration_days": mean_dur,
            })
            print(f"{ds_name} t>{t}d {year}: mean={mean_dur:.2f} d  (n={npx})")

df = pd.DataFrame(rows)
df.to_csv(ROOT / "duration_threshold_mean_results.csv", index=False)

# pivoted year x threshold tables (one per dataset)
for ds_name in DATASETS:
    sub = df[df["dataset"] == ds_name]
    pv = sub.pivot(index="year", columns="threshold_days",
                   values="mean_duration_days").round(2)
    pv.columns = [f">{c}d" for c in pv.columns]
    pv.to_csv(ROOT / f"mean_duration_by_year_{ds_name}.csv")

# trend summary
summary_rows = []
for (ds_name, t), grp in df.groupby(["dataset", "threshold_days"]):
    grp = grp.sort_values("year")
    yrs  = grp["year"].to_numpy(dtype=float)
    vals = grp["mean_duration_days"].to_numpy(dtype=float)

    mk_res = mk.original_test(vals)
    sen    = mk.sens_slope(vals)
    lr     = stats.linregress(yrs, vals)

    summary_rows.append({
        "dataset": ds_name,
        "threshold_days": t,
        "n_years": len(vals),
        "first_year_mean_d": vals[0],
        "last_year_mean_d":  vals[-1],
        "delta_d":           vals[-1] - vals[0],
        "pct_change":        (vals[-1] - vals[0]) / vals[0] * 100.0 if vals[0] else np.nan,
        "mean_d":            vals.mean(),
        "MK_trend":          mk_res.trend,
        "MK_p":              mk_res.p,
        "MK_z":              mk_res.z,
        "Sens_slope_d_per_yr": sen.slope,
        "OLS_slope_d_per_yr":  lr.slope,
        "OLS_intercept":       lr.intercept,
        "OLS_r2":              lr.rvalue ** 2,
        "OLS_p":               lr.pvalue,
    })

summary = pd.DataFrame(summary_rows)
summary.to_csv(ROOT / "duration_threshold_mean_summary.csv", index=False)

with pd.option_context("display.max_columns", None, "display.width", 200):
    print("\n=== MEAN DURATION SUMMARY ===")
    print(summary.to_string(index=False))
