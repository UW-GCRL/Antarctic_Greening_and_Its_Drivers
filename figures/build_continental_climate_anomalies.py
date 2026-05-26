"""Compute annual continental area-weighted means, anomalies, and a per-year
spatial-spread band for the eight climate variables used in the Shapley
analysis.

For each year y, area-weighted across pixels with weights w_p = land_area(p):
    yearly_mean(y)      = sum_p( w_p * var(p,y) ) / sum_p( w_p )
    anomaly(y)          = yearly_mean(y) - mean over years of yearly_mean

For the per-year spatial-spread band, we compute the area-weighted standard
deviation of pixel values around that year's continental mean:
    mean_spatial_std(y) = sqrt( sum_p( w_p * ( var(p,y) - yearly_mean(y) )^2 )
                                / sum_p( w_p ) )
    ci_lo / ci_hi       = anomaly +/- mean_spatial_std

mean_spatial_std measures how heterogeneous the variable is across pixels
within that year (the area-weighted spatial std, the textbook within-year
spatial variability). Useful as a shaded +/- 1 SD envelope on the sparkline.

Saves one CSV per variable in figure1data/.
"""

from pathlib import Path
import numpy as np
import pandas as pd

DATA = Path(r"G:\Hangkai\Anttarctic Vegetation Dynamic\RF_trainning_data\cleared_data")
OUT  = Path(r"G:\Hangkai\Antarctica_Mapping_Data\figure1data")
OUT.mkdir(parents=True, exist_ok=True)

# Column name in raw cleared_all_data_<year>.csv files. Note: "10m_wind_speed"
# in raw is the same variable as "X10m_wind_speed" in the decomposed files
# (R/Python identifier safety prefix).
VARS = [
    "uv_radiation",
    "volumetric_soil_water",
    "temperature_2m",
    "precipitation",
    "solar_radiation",
    "runoff",
    "10m_wind_speed",
    "vapor_pressure_deficit",
]
YEARS = list(range(2002, 2024))

# Short labels for output filenames
SHORT = {
    "uv_radiation": "UV",
    "volumetric_soil_water": "VSW",
    "temperature_2m": "T_2m",
    "precipitation": "Precip",
    "solar_radiation": "Solar",
    "runoff": "Runoff",
    "10m_wind_speed": "Wind",
    "vapor_pressure_deficit": "VPD",
}

# Native unit for each variable's values (yearly_mean, anomaly, std all share
# the same unit since they are linear in the variable). Inferred from ERA5
# documentation cross-checked against actual data magnitudes; the time
# integration period for UV/Solar/Runoff/Precip depends on the upstream
# aggregation in 1-preprocess_rf_data.R and is reported here only as the
# instantaneous/accumulated SI unit.
UNIT = {
    "uv_radiation":            "J m^-2",      # ERA5 downward UV at surface (accumulated)
    "volumetric_soil_water":   "m^3 m^-3",    # fraction 0-1
    "temperature_2m":          "K",           # Kelvin
    "precipitation":           "mm",          # millimetres of water (accumulated)
    "solar_radiation":         "J m^-2",      # ERA5 SSRD (accumulated)
    "runoff":                  "m",           # metres of water (accumulated)
    "10m_wind_speed":          "m s^-1",      # metres per second
    "vapor_pressure_deficit":  "kPa",         # kilopascal
}


def yearly_aw_stats(var):
    """Returns DataFrame with columns:
        year, yearly_mean, anomaly, mean_spatial_std, ci_lo, ci_hi
    """
    by_year = {}
    pix_count = None
    land_area_ref = None
    for y in YEARS:
        f = DATA / f"cleared_all_data_{y}.csv"
        df = pd.read_csv(f, usecols=[var, "land_area"])
        by_year[y] = df[var].to_numpy(dtype=np.float64)
        if pix_count is None:
            pix_count = len(df)
            land_area_ref = df["land_area"].to_numpy(dtype=np.float64)
        else:
            assert len(df) == pix_count, f"row count mismatch in year {y}"

    # Stack as (n_years, n_pixels)
    arr = np.vstack([by_year[y] for y in YEARS])
    w = land_area_ref
    W = w.sum()

    # Per-year area-weighted continental mean
    yearly_mean = (arr * w).sum(axis=1) / W

    # Per-pixel deviation from that year's continental mean (broadcast)
    pixel_dev = arr - yearly_mean[:, None]

    # Area-weighted spatial std of those deviations, per year
    mean_spatial_std = np.sqrt((pixel_dev ** 2 * w).sum(axis=1) / W)

    period_mean = yearly_mean.mean()
    anomaly = yearly_mean - period_mean

    out = pd.DataFrame({
        "year": YEARS,
        "yearly_mean": yearly_mean,
        "anomaly": anomaly,
        "mean_spatial_std": mean_spatial_std,
        "ci_lo": anomaly - mean_spatial_std,
        "ci_hi": anomaly + mean_spatial_std,
    })
    return out


for v in VARS:
    print(f"[compute] {v}  (unit: {UNIT[v]})")
    out = yearly_aw_stats(v)
    out["variable"] = v
    out["unit"] = UNIT[v]
    csv_path = OUT / f"continental_yearly_{SHORT[v]}.csv"
    out[["year", "variable", "unit", "yearly_mean", "anomaly",
         "mean_spatial_std", "ci_lo", "ci_hi"]].to_csv(csv_path, index=False)
    print(f"  saved: {csv_path}")
    print(out[["year", "yearly_mean", "anomaly", "mean_spatial_std",
               "ci_lo", "ci_hi"]].to_string(index=False))
    print()
