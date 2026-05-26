"""Plot the per-region significance matrix across the 256 detrending
counterfactuals -- the figure previously known as main-text Fig. 4 (now
Fig. S23 in the revised manuscript).

Each region is one row block with three sub-panels:
  - Top:    scatter of mean MK p-value across scenarios, log y, colored by
            significance (sig: dark blue, non-sig: orange), 0.05 line.
  - Middle: 8 (variable) x 256 (scenario) dot matrix. A filled dot at
            (v, s) means variable v is detrended in scenario s. Scenarios
            are ordered by mean MK p-value (ascending), so the leftmost
            scenarios are the ones whose green-area trend is most
            statistically significant after detrending.
  - Right:  horizontal bar -- for each variable, fraction of NON-significant
            scenarios (p >= 0.05) that include that variable. Variables
            that appear often in non-significant scenarios are the ones
            whose detrending most consistently removes the trend.

Inputs (computed by compute_mk_all_regions.py):
  Continental_mk_perpixel.csv
  East_mk_perpixel.csv
  West_mk_perpixel.csv
  Peninsula_mk_perpixel.csv

Outputs:
  significance_matrix_perpixel.png / .pdf
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

NEW = Path(r"G:\Hangkai\Antarctica_Mapping_Data\New\DetrendCurves_PerPixelCentered")
OUT_PNG = NEW / "significance_matrix_perpixel.png"
OUT_PDF = NEW / "significance_matrix_perpixel.pdf"

# Variable order in the dot matrix (top -> bottom). Aligned with manuscript
# convention; the variable display labels are short for readability.
VARS = [
    "volumetric_soil_water",
    "precipitation",
    "uv_radiation",
    "temperature_2m",
    "X10m_wind_speed",
    "solar_radiation",
    "runoff",
    "vapor_pressure_deficit",
]
VAR_LABEL = {
    "volumetric_soil_water":  "Soil water",
    "precipitation":          "Precipitation",
    "uv_radiation":           "UV radiation",
    "temperature_2m":         "Air temp. (2 m)",
    "X10m_wind_speed":        "Wind speed (10 m)",
    "10m_wind_speed":         "Wind speed (10 m)",
    "solar_radiation":        "Solar radiation",
    "runoff":                 "Runoff",
    "vapor_pressure_deficit": "VPD",
}

REGIONS = [
    ("Continental", "A. Whole Antarctica"),
    ("East",        "B. East Antarctica"),
    ("West",        "C. West Antarctica"),
    ("Peninsula",   "D. Antarctic Peninsula"),
]

SIG_THRESH = 0.05
SIG_COLOR   = "#2C5E8F"   # dark blue: trend remains significant
NSIG_COLOR  = "#E08A2A"   # orange: trend non-significant after detrending


plt.rcParams.update({
    "font.family":      "DejaVu Sans",
    "font.size":        9.5,
    "axes.titlesize":   10.5,
    "axes.labelsize":   9.5,
    "axes.linewidth":   0.7,
    "xtick.labelsize":  8.5,
    "ytick.labelsize":  8.5,
    "legend.fontsize":  9,
    "axes.spines.top":     False,
    "axes.spines.right":   False,
    "axes.titleweight":   "bold",
    "axes.titlelocation": "left",
    "axes.titlepad":      4,
})


def decode_scenario(sname):
    """Return the set of variables detrended in this scenario name."""
    if sname == "baseline":
        return set()
    return set(p.strip() for p in sname.split(",") if p.strip())


# ---------------------------------------------------------------------------
# Build figure
# ---------------------------------------------------------------------------
n_regions = len(REGIONS)

fig = plt.figure(figsize=(14.0, 3.4 * n_regions + 0.4))
gs_outer = fig.add_gridspec(
    nrows=n_regions, ncols=1,
    hspace=0.55,
    left=0.06, right=0.985, top=0.965, bottom=0.07,
)

for ri, (region, title) in enumerate(REGIONS):
    df = pd.read_csv(NEW / f"{region}_mk_perpixel.csv")
    df = df.sort_values("mean_p_value").reset_index(drop=True)
    n_scen = len(df)
    df["sig"] = df["mean_p_value"] < SIG_THRESH

    # Decode which variables are in each scenario
    has_var = np.zeros((len(VARS), n_scen), dtype=bool)
    for col_idx, (_, row) in enumerate(df.iterrows()):
        in_scen = decode_scenario(row["scenario"])
        for vi, v in enumerate(VARS):
            has_var[vi, col_idx] = (v in in_scen)

    # Frequency of each variable among non-significant scenarios
    nsig_mask = (df["mean_p_value"] >= SIG_THRESH).to_numpy()
    n_nsig = int(nsig_mask.sum())
    freq_nsig = (has_var[:, nsig_mask].sum(axis=1) / max(n_nsig, 1)) * 100.0

    # Inner gridspec for this region's row: 3 sub-axes
    gs_in = gs_outer[ri].subgridspec(
        nrows=2, ncols=2,
        height_ratios=[1.0, 1.6],
        width_ratios=[6.0, 1.0],
        hspace=0.06, wspace=0.06,
    )
    ax_scatter = fig.add_subplot(gs_in[0, 0])
    ax_matrix  = fig.add_subplot(gs_in[1, 0], sharex=ax_scatter)
    ax_freq    = fig.add_subplot(gs_in[1, 1], sharey=ax_matrix)
    # Top-right cell stays empty (no axes)

    # ---- Top: scatter of p-values ----
    xs = np.arange(n_scen)
    colors = [SIG_COLOR if s else NSIG_COLOR for s in df["sig"]]
    ax_scatter.scatter(xs, df["mean_p_value"], c=colors, s=12,
                       edgecolor="none", zorder=3)
    ax_scatter.axhline(SIG_THRESH, color="#666666", ls="--", lw=0.9, zorder=2)
    ax_scatter.set_yscale("log")
    ax_scatter.set_ylabel("MK $p$")
    ax_scatter.tick_params(axis="x", labelbottom=False)
    ax_scatter.set_title(title, fontsize=11)
    ax_scatter.grid(True, which="both", color="0.92", linewidth=0.4, zorder=0)
    ax_scatter.set_axisbelow(True)
    ax_scatter.set_xlim(-1, n_scen)

    # ---- Middle: dot matrix ----
    rows, cols = np.where(has_var)
    # Color each dot by the scenario's significance
    dot_colors = np.where(df.iloc[cols]["sig"].to_numpy(), SIG_COLOR, NSIG_COLOR)
    ax_matrix.scatter(cols, rows, c=dot_colors, s=8,
                      marker="o", edgecolor="none", zorder=3)
    ax_matrix.set_yticks(np.arange(len(VARS)))
    ax_matrix.set_yticklabels([VAR_LABEL[v] for v in VARS])
    ax_matrix.invert_yaxis()
    ax_matrix.set_xlim(-1, n_scen)
    ax_matrix.set_ylim(len(VARS) - 0.5, -0.5)
    ax_matrix.set_xlabel("Scenario  (sorted by MK $p$, ascending)")
    ax_matrix.grid(True, axis="y", color="0.93", linewidth=0.4, zorder=0)
    ax_matrix.set_axisbelow(True)
    ax_matrix.tick_params(axis="x", labelsize=7.5)

    # ---- Right: horizontal bar of variable frequency in non-sig scenarios ----
    ax_freq.barh(
        np.arange(len(VARS)), freq_nsig,
        color=NSIG_COLOR, edgecolor="black", linewidth=0.5, height=0.65,
    )
    for vi, f in enumerate(freq_nsig):
        ax_freq.text(f + 1.0, vi, f"{f:.0f}%", va="center", ha="left",
                     fontsize=7.5, color="0.15")
    ax_freq.set_xlim(0, max(105, freq_nsig.max() + 12))
    ax_freq.invert_yaxis()
    ax_freq.set_xlabel("Frequency in\nnon-sig. scenarios (%)", fontsize=8.5)
    ax_freq.tick_params(axis="y", labelleft=False)
    ax_freq.tick_params(axis="x", labelsize=7.5)
    ax_freq.grid(True, axis="x", color="0.92", linewidth=0.4, zorder=0)
    ax_freq.set_axisbelow(True)

# Shared legend at the bottom
from matplotlib.lines import Line2D
handles = [
    Line2D([0], [0], marker="o", linestyle="", markerfacecolor=SIG_COLOR,
           markeredgecolor=SIG_COLOR, markersize=8,
           label=r"Trend significant after detrending ($p < 0.05$)"),
    Line2D([0], [0], marker="o", linestyle="", markerfacecolor=NSIG_COLOR,
           markeredgecolor=NSIG_COLOR, markersize=8,
           label=r"Trend non-significant ($p \geq 0.05$)"),
]
fig.legend(
    handles=handles,
    loc="lower center", bbox_to_anchor=(0.5, 0.0),
    ncol=2, frameon=False, fontsize=9.5,
)

fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
fig.savefig(OUT_PDF, bbox_inches="tight")
print(f"saved: {OUT_PNG}")
print(f"saved: {OUT_PDF}")
