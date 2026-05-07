"""Combined 3-panel attribution figure: variable-trend table + Shapley
absolute + Shapley percentage."""

import textwrap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import rcParams

from paths import OUTPUT_DIR

REGIONS = ["Continental", "East", "West", "Peninsula"]
VARS    = ["volumetric_soil_water", "precipitation", "uv_radiation",
           "temperature_2m", "X10m_wind_speed", "solar_radiation",
           "runoff", "vapor_pressure_deficit"]
SHORT   = {"volumetric_soil_water": "VSW", "precipitation": "Precip",
           "uv_radiation": "UV", "temperature_2m": "T_2m",
           "X10m_wind_speed": "Wind", "solar_radiation": "Solar",
           "runoff": "Runoff", "vapor_pressure_deficit": "VPD"}
UNITS = {"temperature_2m": "K yr^-1", "uv_radiation": "J m^-2 yr^-1",
         "solar_radiation": "J m^-2 yr^-1",
         "volumetric_soil_water": "m^3 m^-3 yr^-1", "runoff": "m yr^-1",
         "X10m_wind_speed": "m s^-1 yr^-1", "precipitation": "mm yr^-1",
         "vapor_pressure_deficit": "kPa yr^-1"}
COLS = [SHORT[v] for v in VARS]

rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.titleweight": "bold",
    "axes.titlesize": 14,
})

HEADER_BG = "#1f3a5f"
LABEL_BG  = "#1f3a5f"


def soft(color, blend=0.55):
    base = np.array(color)
    base[:3] = base[:3] * (1 - blend) + blend
    return tuple(base)


def fmt_slope_smart(m, s, p):
    am = abs(m)
    if am == 0:
        m_str, s_str = "0", f"{s:.2g}"
    elif am >= 1:
        m_str, s_str = f"{m:+.2f}", f"{s:.2f}"
    elif am >= 1e-2:
        m_str, s_str = f"{m:+.4f}", f"{s:.3f}"
    else:
        m_str, s_str = f"{m:+.2e}", f"{s:.1e}"
    if not np.isnan(p):
        sig = "**" if p < 0.01 else ("*" if p < 0.05 else "")
    else:
        sig = ""
    return f"{m_str}{sig}", f"+/-{s_str}"


def style_cells(table, header_height=0.13, body_height=0.18, header_size=11,
                body_size=11.5, bold_mask=None):
    table.auto_set_font_size(False)
    table.scale(1.0, 1.0)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("white"); cell.set_linewidth(2.5)
        if r == 0:
            cell.set_facecolor(HEADER_BG)
            cell.set_text_props(color="white", weight="bold", size=header_size)
            cell.set_height(header_height)
        elif c == -1:
            cell.set_facecolor(LABEL_BG)
            cell.set_text_props(color="white", weight="bold", size=11)
            cell.set_height(body_height)
        else:
            cell.set_height(body_height)
            cell.set_text_props(size=body_size)
            if bold_mask is not None and bold_mask[r - 1][c]:
                cell.set_text_props(weight="bold", size=body_size)


def render_trend_table(ax, trends_df):
    cmap = plt.cm.RdBu_r
    norms = {}
    for v in VARS:
        vals = trends_df[trends_df.variable == v]["mean_slope_per_pixel"].values
        vmax = max(abs(vals.min()), abs(vals.max()), 1e-12)
        norms[v] = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    cell_text, cell_colors = [], []
    for region in REGIONS:
        row_t, row_c = [], []
        for v in VARS:
            sub = trends_df[(trends_df.region == region) & (trends_df.variable == v)].iloc[0]
            mean_str, std_str = fmt_slope_smart(sub.mean_slope_per_pixel,
                                                sub.std_slope_per_pixel,
                                                sub.mk_p_year_mean)
            row_t.append(f"{mean_str}\n{std_str}")
            row_c.append(soft(cmap(norms[v](sub.mean_slope_per_pixel)), blend=0.5))
        cell_text.append(row_t)
        cell_colors.append(row_c)

    col_labels = [f"{SHORT[v]}\n({UNITS[v]})" for v in VARS]
    table = ax.table(cellText=cell_text, cellColours=cell_colors,
                     rowLabels=REGIONS, colLabels=col_labels,
                     loc="center", cellLoc="center")
    style_cells(table, header_height=0.20, body_height=0.20, header_size=11, body_size=11)
    ax.set_title("(A)  Per-pixel linear trend of climate variables, 2002-2023  "
                 "- area-weighted mean +/- std across pixels  (** p<0.01, * p<0.05 on annual mean)",
                 fontsize=14, pad=14, color="#1f3a5f", loc="left")
    ax.axis("off")


def render_shapley_table(ax, df, fmt, vmax, include_total, title):
    cmap = plt.cm.RdBu_r
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    cell_text, cell_colors, bold_mask = [], [], []
    for region in REGIONS:
        row_t, row_c, row_b = [], [], []
        means = [df.loc[region, f"{c}_mean"] for c in COLS]
        leader = int(np.argmax(np.abs(means)))
        for i, c in enumerate(COLS):
            m = df.loc[region, f"{c}_mean"]
            s = df.loc[region, f"{c}_std"]
            row_t.append(f"{fmt.format(m)}\n+/-{abs(s):.2f}")
            row_c.append(soft(cmap(norm(m)), blend=0.5))
            row_b.append(i == leader)
        if include_total:
            tm = df.loc[region, "TOTAL_mean"]; ts = df.loc[region, "TOTAL_std"]
            row_t.append(f"{fmt.format(tm)}\n+/-{abs(ts):.2f}")
            row_c.append((0.92, 0.94, 0.97, 1.0))
            row_b.append(True)
        cell_text.append(row_t); cell_colors.append(row_c); bold_mask.append(row_b)

    col_labels = COLS + (["TOTAL"] if include_total else [])
    table = ax.table(cellText=cell_text, cellColours=cell_colors,
                     rowLabels=REGIONS, colLabels=col_labels,
                     loc="center", cellLoc="center")
    style_cells(table, header_height=0.16, body_height=0.20,
                header_size=11.5, body_size=11.5, bold_mask=bold_mask)
    ax.set_title(title, fontsize=14, pad=14, color="#1f3a5f", loc="left")
    ax.axis("off")


def main():
    trends_df = pd.read_csv(OUTPUT_DIR / "variable_trends_per_region.csv")
    abs_df    = pd.read_csv(OUTPUT_DIR / "shapley_table_abs_with_std.csv").set_index("region").loc[REGIONS]
    pct_df    = pd.read_csv(OUTPUT_DIR / "shapley_table_pct_with_std.csv").set_index("region").loc[REGIONS]

    FIG_W, FIG_H = 22, 22
    fig = plt.figure(figsize=(FIG_W, FIG_H))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 0.95], hspace=0.50,
                          top=0.97, bottom=0.34, left=0.06, right=0.99)
    ax_a = fig.add_subplot(gs[0]); ax_b = fig.add_subplot(gs[1]); ax_c = fig.add_subplot(gs[2])

    render_trend_table(ax_a, trends_df)
    render_shapley_table(ax_b, abs_df, fmt="{:+.2f}", vmax=15.0, include_total=True,
                         title="(B)  Shapley attribution of greening trend  (km^2 yr^-1, mean +/- std across 5 folds)")
    render_shapley_table(ax_c, pct_df, fmt="{:+.1f}", vmax=55.0, include_total=False,
                         title="(C)  Shapley attribution as % of region's climate-explained trend")

    caption_blocks = [
        ("Panel A - variable trends.",
         "Area-weighted mean (+/- weighted std) of per-pixel OLS slopes for each climate variable in each region, in native units per year. "
         "Stars: Mann-Kendall significance of the area-weighted annual mean (** p<0.01, * p<0.05). "
         "Color: blue = positive trend (variable increasing); red = negative trend; saturation per column."),
        ("Panels B & C - Shapley attribution of the modeled greening trend.",
         "Computed by detrending each variable per pixel (centered on the period mean year) and re-predicting with the random forest; "
         "values across all 256 detrending coalitions are combined via the Shapley formula. "
         "Bolded cells highlight the leading single-variable driver per region."),
        ("Sign convention (Panels B & C).",
         "POSITIVE (blue) = the variable's 2002-2023 trend has, on net, HELPED the greening trend (its removal would slow greening by that amount). "
         "NEGATIVE (red) = the variable's trend has, on net, SLOWED the greening trend (without it, greening would have been faster). "
         "The sign is NOT the variable's own direction of change - it is the product of (variable's trend direction) x (variable's marginal effect on vegetation). "
         "Example: UV declined over 2002-2023 and harms vegetation (negative effect); the two negatives multiply to a positive Shapley value, "
         "meaning UV decline helped greening."),
        ("TOTAL column (Panel B).",
         "= baseline trend - all-8-detrended trend = the climate-explained portion of each region's greening rate."),
    ]

    caption_top = 0.33
    caption_bottom = 0.04
    fig.add_artist(plt.Rectangle((0.025, caption_bottom), 0.965, caption_top - caption_bottom,
                                  transform=fig.transFigure, facecolor="#f7f9fc",
                                  edgecolor="#d0d7e2", linewidth=1.2, zorder=-1))

    WRAP_WIDTH = 185
    y_cursor = caption_top - 0.014
    LINE_H   = 0.0118
    BLOCK_GAP = 0.012
    for title, body in caption_blocks:
        fig.text(0.035, y_cursor, "|", fontsize=14, color="#1f3a5f", weight="bold")
        fig.text(0.045, y_cursor, title, fontsize=12.5, weight="bold", color="#1f3a5f")
        y_cursor -= LINE_H * 1.7
        wrapped = textwrap.wrap(body, width=WRAP_WIDTH)
        for line in wrapped:
            fig.text(0.045, y_cursor, line, fontsize=11, color="#222222")
            y_cursor -= LINE_H
        y_cursor -= BLOCK_GAP

    out = OUTPUT_DIR / "Combined_attribution_figure.png"
    plt.savefig(out, dpi=240, facecolor="white")
    plt.savefig(out.with_suffix(".pdf"), facecolor="white")
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
