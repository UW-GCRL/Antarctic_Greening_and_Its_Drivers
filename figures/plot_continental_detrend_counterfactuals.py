"""Quick continental-scale analysis of the per-pixel-centered detrending results.
Reports:
 - Baseline yearly area + slope
 - Slope of every single-variable detrend scenario (Δ vs baseline)
 - Slope of the all-8-detrended scenario
 - Side-by-side comparison vs the old anchor-at-2002 results (DetrendCurves/)
 - Saves a 5-scenario plot (baseline, detrend T, UV, VSW, UV+VSW) for visual.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

NEW = Path(r"G:\Hangkai\Antarctica_Mapping_Data\New\DetrendCurves_PerPixelCentered")
OLD = Path(r"G:\Hangkai\Antarctica_Mapping_Data\New\DetrendCurves")

VARS = ["temperature_2m", "uv_radiation", "solar_radiation",
        "volumetric_soil_water", "runoff", "X10m_wind_speed",
        "precipitation", "vapor_pressure_deficit"]
ALL_NAME = ",".join(sorted(VARS))


def slope(years, vals):
    a, _ = np.polyfit(years.astype(float), vals, 1)
    return a


df_new = pd.read_csv(NEW / "Continental_detrend_curves.csv")
df_old = pd.read_csv(OLD / "Continental_detrend_curves.csv") if (OLD / "Continental_detrend_curves.csv").exists() else None

print("=" * 80)
print("CONTINENTAL - per-pixel centered detrending - raw analysis")
print("=" * 80)

# Baseline curve
b = df_new[df_new.scenario == "baseline"].sort_values("year")
b_slope = slope(b.year.values, b.mean_area.values)
print(f"\nBaseline yearly area (mean +/- std across 5 folds):")
print(f"{'year':<6} {'mean_area (km²)':>18} {'std':>12}")
for _, row in b.iterrows():
    print(f"{int(row.year):<6} {row.mean_area:>18.2f} {row.std_area:>12.2f}")
print(f"\nBaseline OLS slope = {b_slope:.2f} km2/yr")

# All-detrended
a = df_new[df_new.scenario == ALL_NAME].sort_values("year")
a_slope = slope(a.year.values, a.mean_area.values)
print(f"\nAll-8-detrended OLS slope = {a_slope:.2f} km2/yr  (delta vs baseline = {a_slope - b_slope:+.2f})")

# Single-variable detrend slopes
print(f"\nSingle-variable detrend slopes (sorted by delta from baseline):")
print(f"{'variable':<25} {'slope':>10} {'delta vs baseline':>20}")
rows = []
for v in VARS:
    sub = df_new[df_new.scenario == v].sort_values("year")
    sl = slope(sub.year.values, sub.mean_area.values)
    rows.append((v, sl, sl - b_slope))
rows.sort(key=lambda r: r[2])
for v, sl, dsl in rows:
    print(f"{v:<25} {sl:>10.2f} {dsl:>20.2f}")

# Compare with old method
if df_old is not None:
    print()
    print("=" * 80)
    print("Side-by-side: NEW (per-pixel centered) vs OLD (anchor-at-2002)")
    print("=" * 80)
    print(f"{'scenario':<28} {'NEW slope':>11} {'OLD slope':>11} {'delta (new-old)':>17}")
    for sname in ["baseline"] + VARS + [ALL_NAME]:
        new_sub = df_new[df_new.scenario == sname].sort_values("year")
        old_sub = df_old[df_old.scenario == sname].sort_values("year")
        if len(new_sub) and len(old_sub):
            ns = slope(new_sub.year.values, new_sub.mean_area.values)
            os = slope(old_sub.year.values, old_sub.mean_area.values)
            label = "all-8" if sname == ALL_NAME else sname
            print(f"{label:<28} {ns:>11.2f} {os:>11.2f} {ns - os:>17.2f}")

# Plot
fig, ax = plt.subplots(figsize=(8, 5.2))
plot_specs = [
    ("baseline",                              "Baseline (no detrending)",   "#000000", "-"),
    ("temperature_2m",                        "Detrend Temperature",        "#d62728", "-"),
    ("uv_radiation",                          "Detrend UV radiation",       "#ff7f0e", "-"),
    ("volumetric_soil_water",                 "Detrend Soil moisture",      "#1f77b4", "-"),
    ("uv_radiation,volumetric_soil_water",    "Detrend UV + Soil moisture", "#9467bd", "--"),
    (ALL_NAME,                                "Detrend ALL 8 variables",    "#2ca02c", ":"),
]
for sname, label, color, ls in plot_specs:
    sub = df_new[df_new.scenario == sname].sort_values("year")
    if not len(sub):
        continue
    yrs, m, s = sub.year.values, sub.mean_area.values, sub.std_area.values
    ax.plot(yrs, m, color=color, ls=ls, lw=1.8, label=label)
    ax.fill_between(yrs, m - s, m + s, color=color, alpha=0.18, lw=0)

ax.set_xlabel("Year")
ax.set_ylabel("Continental green area (km²)")
ax.legend(loc="upper left", frameon=False, fontsize=9)
ax.grid(True, alpha=0.3)
plt.tight_layout()
out = NEW / "Continental_detrend_perpixel.png"
plt.savefig(out, dpi=200)
plt.savefig(out.with_suffix(".pdf"))
print(f"\nSaved plot: {out}")
