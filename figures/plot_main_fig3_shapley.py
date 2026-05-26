"""Figure 3 -- Shapley attribution figures.

Produces two variants:
  1. figure3_regional_shapley.{png,pdf}                 -- 1 x 3 (East / West / Peninsula)
  2. figure3_continental_and_regional_shapley.{png,pdf} -- 2 x 2 (Continental + 3 regions)

For each panel: bars are sorted by phi ascending; positive bars (the
variable's 2002-2023 trend helped greening) are green; negative bars
(slowed greening) are red. Error bars are 1 std across 5 RF folds.
Sigma phi (climate-explained portion) is annotated in each panel.
"""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

DATA = Path(r"G:\Hangkai\Antarctica_Mapping_Data\figure1data")
OUT_PNG = DATA / "figure3_regional_shapley.png"
OUT_PDF = DATA / "figure3_regional_shapley.pdf"
OUT_PNG_4 = DATA / "figure3_continental_and_regional_shapley.png"
OUT_PDF_4 = DATA / "figure3_continental_and_regional_shapley.pdf"

DISPLAY = {
    "volumetric_soil_water":   r"Soil water",
    "precipitation":           r"Precipitation",
    "uv_radiation":            r"UV radiation",
    "temperature_2m":          r"Air temp. (2 m)",
    "X10m_wind_speed":         r"Wind speed (10 m)",
    "10m_wind_speed":          r"Wind speed (10 m)",
    "solar_radiation":         r"Solar radiation",
    "runoff":                  r"Runoff",
    "vapor_pressure_deficit":  r"VPD",
}

REGIONS_3 = [
    ("East",      "East Antarctica",     DATA / "East_shapley.csv"),
    ("West",      "West Antarctica",     DATA / "West_shapley.csv"),
    ("Peninsula", "Antarctic Peninsula", DATA / "Peninsula_shapley.csv"),
]

REGIONS_4 = [
    ("Continental", "Continental",         DATA / "Continental_shapley.csv"),
    ("East",        "East Antarctica",     DATA / "East_shapley.csv"),
    ("West",        "West Antarctica",     DATA / "West_shapley.csv"),
    ("Peninsula",   "Antarctic Peninsula", DATA / "Peninsula_shapley.csv"),
]

POS_COLOR = "#2E8B57"
NEG_COLOR = "#C8493B"


def load_panels(region_list):
    panels = []
    for key, title, fp in region_list:
        df = pd.read_csv(fp)
        df = df.sort_values("shapley_mean", ascending=True).reset_index(drop=True)
        panels.append((key, title, df))
    return panels


def shared_xlim(panels):
    all_means = pd.concat([d[["shapley_mean", "shapley_std"]] for _, _, d in panels])
    xmin = (all_means["shapley_mean"] - all_means["shapley_std"]).min()
    xmax = (all_means["shapley_mean"] + all_means["shapley_std"]).max()
    pad  = 0.10 * (xmax - xmin)
    return (xmin - pad, xmax + pad)


def draw_panel(ax, df, title, panel_label, xlim, label_fs=10, title_fs=12):
    y = range(len(df))
    means = df["shapley_mean"].to_numpy()
    stds  = df["shapley_std"].to_numpy()
    colors = [POS_COLOR if m >= 0 else NEG_COLOR for m in means]
    labels = [DISPLAY.get(v, v) for v in df["variable"]]

    ax.barh(
        y, means, xerr=stds,
        color=colors, edgecolor="black", linewidth=0.6, height=0.7,
        error_kw=dict(ecolor="black", elinewidth=1.0, capsize=3.0),
        zorder=3,
    )

    ax.axvline(0, color="black", linewidth=0.8, zorder=2)
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=label_fs)
    ax.set_xlim(*xlim)
    ax.set_xlabel(r"Shapley value $\varphi$  (km$^2$ yr$^{-1}$)", fontsize=label_fs)
    ax.set_title(f"{panel_label}. {title}", loc="left",
                 fontsize=title_fs, fontweight="bold")
    ax.tick_params(axis="x", labelsize=label_fs - 1)
    ax.grid(axis="x", color="0.85", linewidth=0.6, zorder=1)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    for yi, (m, s) in enumerate(zip(means, stds)):
        if m < 0:
            x_text = m - s - 0.02 * (xlim[1] - xlim[0])
            ha = "right"
        else:
            x_text = m + s + 0.02 * (xlim[1] - xlim[0])
            ha = "left"
        ax.text(x_text, yi, f"{m:+.2f}", va="center", ha=ha,
                fontsize=label_fs - 1.5, color="0.2")

    total = means.sum()
    ax.text(
        0.97, 0.05,
        rf"$\Sigma\varphi = {total:+.2f}$ km$^2$ yr$^{{-1}}$",
        transform=ax.transAxes,
        ha="right", va="bottom",
        fontsize=label_fs - 0.5, color="0.15",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                  edgecolor="0.6", linewidth=0.6, alpha=0.95),
        zorder=5,
    )


def add_legend(fig):
    legend_handles = [
        Patch(facecolor=POS_COLOR, edgecolor="black",
              label=r"Variable trend accelerated green area expansion ($\varphi > 0$)"),
        Patch(facecolor=NEG_COLOR, edgecolor="black",
              label=r"Variable trend slowed green area expansion ($\varphi < 0$)"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.005),
        ncol=2, frameon=False, fontsize=10,
    )


# =========================================================================
# Variant 1: 1 x 3 (East / West / Peninsula)
# =========================================================================
panels3 = load_panels(REGIONS_3)
xlim3 = shared_xlim(panels3)

fig, axes = plt.subplots(1, 3, figsize=(14.0, 5.0), sharey=False)
plt.subplots_adjust(left=0.085, right=0.985, top=0.92, bottom=0.18, wspace=0.55)

for ax, (_, title, df), pl in zip(axes, panels3, ["A", "B", "C"]):
    draw_panel(ax, df, title, pl, xlim3)

add_legend(fig)
fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
fig.savefig(OUT_PDF, bbox_inches="tight")
print(f"saved: {OUT_PNG}")
print(f"saved: {OUT_PDF}")
plt.close(fig)


# =========================================================================
# Variant 2: 2 x 2 (Continental + East / West / Peninsula)
# =========================================================================
panels4 = load_panels(REGIONS_4)
xlim4 = shared_xlim(panels4)

fig, axes = plt.subplots(2, 2, figsize=(13.0, 9.0), sharey=False)
plt.subplots_adjust(left=0.085, right=0.985, top=0.95,
                    bottom=0.09, wspace=0.55, hspace=0.42)

for ax, (_, title, df), pl in zip(axes.flat, panels4, ["A", "B", "C", "D"]):
    draw_panel(ax, df, title, pl, xlim4)

add_legend(fig)
fig.savefig(OUT_PNG_4, dpi=300, bbox_inches="tight")
fig.savefig(OUT_PDF_4, bbox_inches="tight")
print(f"saved: {OUT_PNG_4}")
print(f"saved: {OUT_PDF_4}")
plt.close(fig)
