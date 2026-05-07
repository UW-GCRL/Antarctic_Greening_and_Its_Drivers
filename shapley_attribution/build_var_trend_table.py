"""Compute area-weighted per-pixel OLS trend for each climate variable in each
region (Continental, East, West, Peninsula), in the variable's native units
per year. Saves a CSV and a publication-quality figure with Mann-Kendall
significance stars on the area-weighted annual mean."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pymannkendall as mk

from paths import DATA, OUTPUT_DIR

VARS = ["temperature_2m", "uv_radiation", "solar_radiation",
        "volumetric_soil_water", "runoff", "X10m_wind_speed",
        "precipitation", "vapor_pressure_deficit"]
SHORT = {"temperature_2m": "T_2m", "uv_radiation": "UV",
         "solar_radiation": "Solar", "volumetric_soil_water": "VSW",
         "runoff": "Runoff", "X10m_wind_speed": "Wind",
         "precipitation": "Precip", "vapor_pressure_deficit": "VPD"}
UNITS = {
    "temperature_2m": "K yr^-1",
    "uv_radiation":   "J m^-2 yr^-1",
    "solar_radiation": "J m^-2 yr^-1",
    "volumetric_soil_water": "m^3 m^-3 yr^-1",
    "runoff": "m yr^-1",
    "X10m_wind_speed": "m s^-1 yr^-1",
    "precipitation": "mm yr^-1",
    "vapor_pressure_deficit": "kPa yr^-1",
}
REGIONS = ["Continental", "East", "West", "Peninsula"]


def load_baseline():
    df = pd.read_csv(DATA / "detrended_temperature_2m.csv")
    df_uv = pd.read_csv(DATA / "detrended_uv_radiation.csv",
                        usecols=["temperature_2m_temporal", "temperature_2m_residual"])
    df["temperature_2m_temporal"] = df_uv["temperature_2m_temporal"].values
    df["temperature_2m_residual"] = df_uv["temperature_2m_residual"].values
    return df


def per_pixel_slopes(full, pixel_idx, sum_yc2, count_p, y_centered):
    n_pixels = sum_yc2.size
    sum_p = np.zeros(n_pixels)
    np.add.at(sum_p, pixel_idx, full)
    mn = sum_p / count_p
    fc = full - mn[pixel_idx]
    numer = np.zeros(n_pixels)
    np.add.at(numer, pixel_idx, y_centered * fc)
    return numer / sum_yc2


def fmt_slope(m, s, p):
    am = abs(m)
    if am == 0 or am >= 1e-2:
        m_str, s_str = f"{m:+.3g}", f"{s:.2g}"
    else:
        m_str, s_str = f"{m:+.2e}", f"{s:.1e}"
    sig = "**" if (not np.isnan(p) and p < 0.01) else ("*" if (not np.isnan(p) and p < 0.05) else "")
    return f"{m_str}{sig}\n+/-{s_str}"


def main():
    df = load_baseline()
    year = df["year"].values
    y_mean = (year.min() + year.max()) / 2.0

    records = []
    for region in REGIONS:
        if region == "Continental":
            mask = np.ones(len(df), dtype=bool)
        else:
            mask = df["Regions"].values == region
        sub = df.loc[mask]
        pixel = sub["pixel"].values
        yr = sub["year"].values
        yc = (yr - y_mean).astype(np.float64)
        land = sub["land_area"].values
        unique, pidx = np.unique(pixel, return_inverse=True)
        n_pix = len(unique)
        sum_yc2 = np.zeros(n_pix); np.add.at(sum_yc2, pidx, yc ** 2)
        count_p = np.zeros(n_pix); np.add.at(count_p, pidx, 1.0)
        sum_land_p = np.zeros(n_pix); np.add.at(sum_land_p, pidx, land)
        pixel_area = sum_land_p / count_p
        sum_area = pixel_area.sum()

        for v in VARS:
            sp = sub[f"{v}_spatial"].values
            tp = sub[f"{v}_temporal"].values
            rs = sub[f"{v}_residual"].values
            full = sp + tp + rs
            slopes = per_pixel_slopes(full, pidx, sum_yc2, count_p, yc)
            mean_slope = float(np.sum(slopes * pixel_area) / sum_area)
            std_slope  = float(np.sqrt(np.sum(pixel_area * (slopes - mean_slope) ** 2) / sum_area))
            df_v = pd.DataFrame({"year": yr, "full": full, "land": land})
            num_y = df_v.groupby("year").apply(lambda g: (g["full"] * g["land"]).sum() / g["land"].sum())
            try:
                p = mk.original_test(num_y.values).p
            except Exception:
                p = np.nan
            records.append({
                "region": region, "variable": v,
                "mean_slope_per_pixel": mean_slope,
                "std_slope_per_pixel":  std_slope,
                "mk_p_year_mean":       p,
            })

    trends = pd.DataFrame(records)
    trends.to_csv(OUTPUT_DIR / "variable_trends_per_region.csv", index=False)
    print(f"Saved: {OUTPUT_DIR / 'variable_trends_per_region.csv'}")

    # Figure
    fig, ax = plt.subplots(figsize=(15, 6.5))
    cmap = plt.cm.RdBu_r
    norm_per_var = {}
    for v in VARS:
        vals = trends[trends["variable"] == v]["mean_slope_per_pixel"].values
        vmax = max(abs(vals.min()), abs(vals.max()), 1e-12)
        norm_per_var[v] = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    cell_text, cell_colors = [], []
    for region in REGIONS:
        row_t, row_c = [], []
        for v in VARS:
            sub = trends[(trends.region == region) & (trends.variable == v)].iloc[0]
            row_t.append(fmt_slope(sub.mean_slope_per_pixel, sub.std_slope_per_pixel, sub.mk_p_year_mean))
            base = np.array(cmap(norm_per_var[v](sub.mean_slope_per_pixel)))
            base[:3] = base[:3] * 0.45 + 0.55
            row_c.append(tuple(base))
        cell_text.append(row_t); cell_colors.append(row_c)

    col_labels = [f"{SHORT[v]}\n({UNITS[v]})" for v in VARS]
    table = ax.table(cellText=cell_text, cellColours=cell_colors,
                     rowLabels=REGIONS, colLabels=col_labels,
                     loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.0, 3.2)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("white"); cell.set_linewidth(2)
        if r == 0:
            cell.set_facecolor("#2c3e50")
            cell.set_text_props(color="white", weight="bold", size=10)
            cell.set_height(0.13)
        if c == -1:
            cell.set_facecolor("#2c3e50")
            cell.set_text_props(color="white", weight="bold", size=10)

    ax.set_title("Per-pixel linear trend of climate variables, 2002-2023  "
                 "(area-weighted mean +/- std across pixels in each region)",
                 fontsize=12, weight="bold", pad=14, color="#222222")
    ax.axis("off")

    caption = (
        "Each cell shows the area-weighted mean per-pixel OLS slope of the reconstructed signal\n"
        "(spatial + temporal + residual) across pixels in the region, in the variable's native units per year.\n"
        "The +/- value is the area-weighted standard deviation of pixel slopes within the region.\n"
        "Stars: Mann-Kendall significance of the region's area-weighted annual-mean time series\n"
        "  (** p < 0.01, * p < 0.05).\n"
        "Color: blue = positive trend; red = negative trend; saturation per column."
    )
    fig.text(0.045, 0.005, caption, fontsize=9.5, va="bottom", linespacing=1.55, color="#333333")

    plt.tight_layout(rect=[0, 0.18, 1, 0.99])
    out = OUTPUT_DIR / "Variable_trends_per_region.png"
    plt.savefig(out, dpi=220, bbox_inches="tight", facecolor="white")
    plt.savefig(out.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
