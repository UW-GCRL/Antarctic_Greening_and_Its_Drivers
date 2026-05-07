"""Render the Shapley tables (absolute and percentage) as a 2-panel
publication-quality figure."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import rcParams

from paths import OUTPUT_DIR

REGIONS = ["Continental", "East", "West", "Peninsula"]
COLS = ["VSW", "Precip", "UV", "T_2m", "Wind", "Solar", "Runoff", "VPD"]

rcParams["font.family"] = "DejaVu Sans"
rcParams["axes.titleweight"] = "bold"


def cell_text(mean, std, fmt):
    return f"{fmt.format(mean)}\n+/- {abs(std):.2f}"


def render_table(ax, df, value_cols, fmt, vmax, title, include_total=True):
    cmap = plt.cm.RdBu_r
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    cell_text_arr, cell_colors_arr, bold_mask = [], [], []
    for region in REGIONS:
        row_text, row_color, row_bold = [], [], []
        means_for_region = [df.loc[region, f"{c}_mean"] for c in value_cols]
        leader_idx = int(np.argmax(np.abs(means_for_region)))
        for i, c in enumerate(value_cols):
            m = df.loc[region, f"{c}_mean"]
            s = df.loc[region, f"{c}_std"]
            row_text.append(cell_text(m, s, fmt))
            base = np.array(cmap(norm(m)))
            base[:3] = base[:3] * 0.4 + 0.6
            row_color.append(tuple(base))
            row_bold.append(i == leader_idx)
        if include_total:
            tm = df.loc[region, "TOTAL_mean"]
            ts = df.loc[region, "TOTAL_std"]
            row_text.append(cell_text(tm, ts, fmt))
            row_color.append((0.92, 0.92, 0.94, 1.0))
            row_bold.append(True)
        cell_text_arr.append(row_text)
        cell_colors_arr.append(row_color)
        bold_mask.append(row_bold)

    col_labels = list(value_cols) + (["TOTAL"] if include_total else [])
    table = ax.table(cellText=cell_text_arr, cellColours=cell_colors_arr,
                     rowLabels=REGIONS, colLabels=col_labels,
                     loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.0, 2.6)

    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("white"); cell.set_linewidth(2)
        if r == 0:
            cell.set_facecolor("#2c3e50")
            cell.set_text_props(color="white", weight="bold", size=10)
            cell.set_height(0.10)
        elif c == -1:
            cell.set_facecolor("#2c3e50")
            cell.set_text_props(color="white", weight="bold", size=10)
        else:
            if bold_mask[r - 1][c]:
                cell.set_text_props(weight="bold", size=10)

    ax.set_title(title, fontsize=12, weight="bold", pad=14, color="#222222")
    ax.axis("off")


def main():
    abs_df = pd.read_csv(OUTPUT_DIR / "shapley_table_abs_with_std.csv").set_index("region").loc[REGIONS]
    pct_df = pd.read_csv(OUTPUT_DIR / "shapley_table_pct_with_std.csv").set_index("region").loc[REGIONS]

    fig = plt.figure(figsize=(14, 9.5))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.05, 1.0], hspace=0.35)
    ax_abs = fig.add_subplot(gs[0])
    ax_pct = fig.add_subplot(gs[1])

    render_table(ax_abs, abs_df, COLS, fmt="{:+.2f}", vmax=15.0, include_total=True,
                 title="Shapley attribution of the 2002-2023 greening trend  (km^2 yr^-1, mean +/- std across 5 folds)")
    render_table(ax_pct, pct_df, COLS, fmt="{:+.1f}", vmax=55.0, include_total=False,
                 title="Shapley attribution as % of each region's climate-explained trend")

    caption = (
        "Sign convention: positive (blue) values mean the variable's 2002-2023 trend has, on net, helped\n"
        "the greening trend (its removal would slow greening by that amount). Negative (red) values mean\n"
        "the variable's trend has, on net, slowed the greening trend. The sign is NOT the direction of\n"
        "change of the variable itself; it is the product of (variable trend direction) x (variable\n"
        "marginal effect on vegetation). Bolded cells: the leading single-variable driver in each region."
    )
    fig.text(0.045, 0.005, caption, fontsize=9.5, va="bottom", linespacing=1.55, color="#333333")

    plt.tight_layout(rect=[0, 0.13, 1, 0.99])
    out = OUTPUT_DIR / "Shapley_attribution_table.png"
    plt.savefig(out, dpi=220, bbox_inches="tight", facecolor="white")
    plt.savefig(out.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
