"""Science-style revision of the duration-threshold robustness figures.

Produces TWO separate 4-panel figures (no figure title; titles go in the
caption):

  duration_threshold_trends.{png,pdf}        -- GREEN AREA
  duration_threshold_mean_trends.{png,pdf}   -- MEAN GREENNESS DURATION

Each 2 x 2 panel:
  A -- annual time series under thresholds >=1..>=14 d, first-last method
  B -- same, daily-count method
  C -- Mann-Kendall p-value vs threshold for both methods (log scale)
  D -- Sen's slope vs threshold for both methods

Reads:
  duration_threshold_trend_results.csv     -- AREA, per-year per-threshold
  duration_threshold_trend_summary.csv     -- AREA, per-threshold MK + Sen
  duration_threshold_mean_results.csv      -- DURATION, per-year per-threshold
  duration_threshold_mean_summary.csv      -- DURATION, per-threshold MK + Sen
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap


def truncate_cmap(cmap, lo=0.0, hi=1.0, n=256):
    """Return a new colormap that uses only the [lo, hi] slice of `cmap`."""
    xs = np.linspace(lo, hi, n)
    return LinearSegmentedColormap.from_list(
        f"{cmap.name}_trunc", cmap(xs), N=n,
    )


# Low-saturation sequential palette: prefer seaborn's "mako" (dark teal to
# pale mint), truncated so the lightest end is still readable on white
# (mako's far-light end is near-white). Fall back to YlGnBu_r if seaborn
# is not installed.
try:
    import seaborn as sns
    _BASE = sns.color_palette("mako", as_cmap=True)
except Exception:
    _BASE = cm.get_cmap("YlGnBu")
# seaborn's mako runs DARK (near-black) -> LIGHT (pale mint) across the
# 0..1 index. We trim BOTH ends: lo=0.20 avoids near-black, hi=0.65 stops
# before the colormap turns pastel so every curve remains readable on white.
SEQ_CMAP = truncate_cmap(_BASE, lo=0.20, hi=0.65)

ROOT = Path(
    r"G:\Hangkai\Anttarctic Vegetation Dynamic\Version_2_data\Antarctica_Vegetation_Duration"
)

DS_LABEL = {
    "first_last":  "First-last day method",
    "daily_count": "Daily-count method",
}
DS_COLOR = {
    "first_last":  "#3B6FB6",   # muted blue
    "daily_count": "#E08A2A",   # muted orange
}
DS_MARKER = {
    "first_last":  "o",
    "daily_count": "s",
}

plt.rcParams.update({
    "font.family":        "DejaVu Sans",
    "font.size":          10,
    "axes.titlesize":     11.5,
    "axes.labelsize":     10,
    "axes.linewidth":     0.8,
    "xtick.labelsize":    9,
    "ytick.labelsize":    9,
    "legend.fontsize":    9,
    "axes.spines.top":     False,
    "axes.spines.right":   False,
    "axes.titleweight":   "bold",
    "axes.titlelocation": "left",
    "axes.titlepad":      6,
})


def make_figure(results, summary, ycol, slope_col, ylabel_ts, ylabel_slope,
                out_png, out_pdf, *, share_top_ylim=True):
    """Build a 2x2 robustness figure for one quantity."""
    thresholds = sorted(results["threshold_days"].unique())
    norm = plt.Normalize(min(thresholds), max(thresholds))
    cmap = SEQ_CMAP

    fig = plt.figure(figsize=(12.0, 8.4))
    gs = fig.add_gridspec(
        nrows=2, ncols=2,
        hspace=0.36, wspace=0.22,
        left=0.075, right=0.895, top=0.96, bottom=0.075,
    )
    ax_A = fig.add_subplot(gs[0, 0])
    ax_B = fig.add_subplot(gs[0, 1])
    ax_C = fig.add_subplot(gs[1, 0])
    ax_D = fig.add_subplot(gs[1, 1])

    # ------- A / B: per-threshold time series ----------
    def plot_ts(ax, ds, panel_label):
        sub = results[results["dataset"] == ds]
        for t in thresholds:
            s = sub[sub["threshold_days"] == t].sort_values("year")
            ax.plot(
                s["year"], s[ycol],
                color=cmap(norm(t)), lw=1.4,
                marker="o", ms=3.0, markeredgewidth=0,
            )
        ax.set_xlabel("Year")
        ax.set_ylabel(ylabel_ts)
        ax.set_title(f"{panel_label}. {DS_LABEL[ds]}")
        ax.grid(axis="y", color="0.88", linewidth=0.5, zorder=0)
        ax.set_axisbelow(True)
        ax.set_xlim(2001.5, 2023.5)

    plot_ts(ax_A, "first_last",  "A")
    plot_ts(ax_B, "daily_count", "B")
    if share_top_ylim:
        ymin = min(ax_A.get_ylim()[0], ax_B.get_ylim()[0])
        ymax = max(ax_A.get_ylim()[1], ax_B.get_ylim()[1])
        ax_A.set_ylim(ymin, ymax)
        ax_B.set_ylim(ymin, ymax)

    # ------- shared colorbar to the right of the top row ----------
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    posB = ax_B.get_position()
    cax = fig.add_axes([
        posB.x1 + 0.012,
        posB.y0,
        0.016,
        posB.y1 - posB.y0,
    ])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label(r"Duration threshold ($\geq T$ days)", fontsize=9.5)
    cbar.set_ticks([1, 2, 4, 6, 8, 10, 12, 14])
    cbar.outline.set_linewidth(0.6)
    cbar.ax.tick_params(labelsize=8.5, width=0.6)

    # ------- C: MK p-value vs threshold ----------
    for ds in ("first_last", "daily_count"):
        s = summary[summary["dataset"] == ds].sort_values("threshold_days")
        ax_C.plot(
            s["threshold_days"], s["MK_p"],
            marker=DS_MARKER[ds], ms=6, lw=1.6,
            color=DS_COLOR[ds], markeredgecolor="black",
            markeredgewidth=0.5, label=DS_LABEL[ds], zorder=3,
        )
    ax_C.axhline(0.05, color="#C8493B", ls="--", lw=1.0,
                 label=r"$p = 0.05$", zorder=2)
    ax_C.set_yscale("log")
    ax_C.set_xlabel(r"Duration threshold ($\geq T$ days)")
    ax_C.set_ylabel(r"Mann-Kendall $p$-value")
    ax_C.set_title("C. Trend significance vs threshold")
    ax_C.grid(True, which="both", color="0.90", linewidth=0.5, zorder=0)
    ax_C.set_axisbelow(True)
    ax_C.set_xticks(thresholds)
    ax_C.legend(loc="best", frameon=False)

    # ------- D: Sen's slope vs threshold ----------
    for ds in ("first_last", "daily_count"):
        s = summary[summary["dataset"] == ds].sort_values("threshold_days")
        ax_D.plot(
            s["threshold_days"], s[slope_col],
            marker=DS_MARKER[ds], ms=6, lw=1.6,
            color=DS_COLOR[ds], markeredgecolor="black",
            markeredgewidth=0.5, label=DS_LABEL[ds], zorder=3,
        )
    ax_D.set_xlabel(r"Duration threshold ($\geq T$ days)")
    ax_D.set_ylabel(ylabel_slope)
    ax_D.set_title("D. Trend magnitude vs threshold")
    ax_D.grid(axis="y", color="0.90", linewidth=0.5, zorder=0)
    ax_D.set_axisbelow(True)
    ax_D.set_xticks(thresholds)
    ax_D.legend(loc="best", frameon=False)

    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    print(f"saved: {out_png}")
    print(f"saved: {out_pdf}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 1 -- Green area robustness
# ---------------------------------------------------------------------------
area_results = pd.read_csv(ROOT / "duration_threshold_trend_results.csv")
area_summary = pd.read_csv(ROOT / "duration_threshold_trend_summary.csv")

make_figure(
    area_results, area_summary,
    ycol="area_km2",
    slope_col="Sens_slope_km2_per_yr",
    ylabel_ts=r"Green area (km$^2$)",
    ylabel_slope=r"Sen's slope (km$^2$ yr$^{-1}$)",
    out_png=ROOT / "duration_threshold_trends.png",
    out_pdf=ROOT / "duration_threshold_trends.pdf",
    share_top_ylim=True,
)


# ---------------------------------------------------------------------------
# Figure 2 -- Mean greenness duration robustness
# ---------------------------------------------------------------------------
dur_results = pd.read_csv(ROOT / "duration_threshold_mean_results.csv")
dur_summary = pd.read_csv(ROOT / "duration_threshold_mean_summary.csv")

# First-last and daily-count duration sit on very different absolute scales
# (~40-90 d vs ~10-40 d), so we do NOT share y-limits between A and B here.
make_figure(
    dur_results, dur_summary,
    ycol="mean_duration_days",
    slope_col="Sens_slope_d_per_yr",
    ylabel_ts="Mean duration (days)",
    ylabel_slope=r"Sen's slope (days yr$^{-1}$)",
    out_png=ROOT / "duration_threshold_mean_trends.png",
    out_pdf=ROOT / "duration_threshold_mean_trends.pdf",
    share_top_ylim=False,
)
