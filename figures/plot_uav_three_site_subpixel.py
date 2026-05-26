"""Science-style revision of three_site_uav_only.png.

Six panels (no figure title; title belongs in the caption):
  A, B, C -- UAV NDVI map with moss extent overlaid (one per site)
  D       -- UAV moss area (m^2) per site
  E       -- UAV moss area as a fraction of one MCD43A4 pixel (%)
  F       -- UAV footprint vs MCD43A4 500-m pixel, drawn to scale

NDVI threshold = 0.2, native UAV resolution = 3 cm. Inputs are the
per-site NPZ files (`*_uav_ndvi.npz`) and the site summary JSON.
"""

from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

ROOT = Path(r"G:\Hangkai\Anttarctic Vegetation Dynamic\AAS_4046")
OUT_PNG = ROOT / "three_site_uav_only.png"
OUT_PDF = ROOT / "three_site_uav_only.pdf"

SITES = ["ASPA135", "RedShed", "RobinsonRidge"]
SITE_DISPLAY = {
    "ASPA135":       "ASPA 135",
    "RedShed":       "Red Shed",
    "RobinsonRidge": "Robinson Ridge",
}
SITE_COLORS = {
    "ASPA135":       "#3B6FB6",   # muted blue
    "RedShed":       "#E08A2A",   # muted orange
    "RobinsonRidge": "#3F8E3F",   # muted green
}
MOSS_COLOR = "#1F6B2E"            # deep green for moss extent

plt.rcParams.update({
    "font.family":      "DejaVu Sans",
    "font.size":        10,
    "axes.titlesize":   11,
    "axes.labelsize":   10,
    "axes.linewidth":   0.8,
    "xtick.labelsize":  9,
    "ytick.labelsize":  9,
    "legend.fontsize":  9,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.titleweight": "bold",
    "axes.titlelocation": "left",
    "axes.titlepad":    6,
    "savefig.dpi":      300,
})


# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
with open(ROOT / "three_site_uav_summary.json") as f:
    summary = json.load(f)

site_data = {}
for site in SITES:
    npz = np.load(ROOT / f"{site}_uav_ndvi.npz", allow_pickle=True)
    site_data[site] = {
        "ndvi":   npz["ndvi"],
        "moss":   npz["moss_mask"].astype(bool),
        "bounds": npz["bounds"],
        "res":    npz["res"],
    }


# ---------------------------------------------------------------------------
# Figure layout
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(13.0, 8.0))
gs = fig.add_gridspec(
    nrows=2, ncols=3,
    height_ratios=[1.65, 1.0],
    hspace=0.30, wspace=0.20,
    left=0.06, right=0.985, top=0.97, bottom=0.085,
)
ax_maps = [fig.add_subplot(gs[0, k]) for k in range(3)]
# Anchor each map axes to the TOP of its grid cell so the panel tops align
# regardless of each site's UAV aspect ratio.
for ax in ax_maps:
    ax.set_anchor("N")
ax_D    = fig.add_subplot(gs[1, 0])
ax_E    = fig.add_subplot(gs[1, 1])
ax_F    = fig.add_subplot(gs[1, 2])


# ---------------------------------------------------------------------------
# A / B / C -- UAV moss-extent maps
# ---------------------------------------------------------------------------
ndvi_cmap = LinearSegmentedColormap.from_list(
    "ndvi_soft", ["#e8e3d1", "#d6cca5", "#bdb978", "#8aa84a", "#3F8E3F"], N=256,
)

panel_labels_top = ["A", "B", "C"]
for ax, site, lbl in zip(ax_maps, SITES, panel_labels_top):
    d = site_data[site]
    info = summary[site]

    # downsample for plotting (3 cm full-res is too big)
    target_max_dim = 900
    H, W = d["ndvi"].shape
    step = max(1, int(np.ceil(max(H, W) / target_max_dim)))
    ndvi_ds = d["ndvi"][::step, ::step]
    moss_ds = d["moss"][::step, ::step]

    xmin, ymin, xmax, ymax = d["bounds"]

    # background NDVI (mask out NaN / out-of-frame)
    base = np.where(np.isnan(ndvi_ds), np.nan, np.clip(ndvi_ds, -0.2, 0.6))
    ax.imshow(
        base, extent=(xmin, xmax, ymin, ymax), origin="upper",
        cmap=ndvi_cmap, vmin=-0.2, vmax=0.6, interpolation="nearest",
    )
    # moss overlay (transparent everywhere else)
    moss_overlay = np.where(moss_ds, 1.0, np.nan)
    ax.imshow(
        moss_overlay, extent=(xmin, xmax, ymin, ymax), origin="upper",
        cmap=plt.matplotlib.colors.ListedColormap([MOSS_COLOR]),
        vmin=0, vmax=1, interpolation="nearest",
    )

    # show ground extent in meters, relative to top-left corner
    xspan = xmax - xmin
    yspan = ymax - ymin
    ax.set_xticks([xmin, xmin + xspan])
    ax.set_xticklabels(["0", f"{xspan:.0f} m"])
    ax.set_yticks([ymin, ymin + yspan])
    ax.set_yticklabels([f"{yspan:.0f} m", "0"])
    ax.set_aspect("equal")

    ax.set_title(
        f"{lbl}. {SITE_DISPLAY[site]}",
        fontsize=12,
    )
    # Stats text is added below in figure coordinates after all axes are
    # placed, so the three blocks share a uniform y position regardless of
    # individual map aspect ratios (otherwise C, the tallest map, would
    # push its stats text into panel F).
    pass

# Shared legend handles for the three UAV maps
legend_handles = [
    Line2D([0], [0], marker="s", linestyle="",
           markerfacecolor=MOSS_COLOR, markeredgecolor=MOSS_COLOR,
           markersize=10, label="Moss extent (UAV NDVI > 0.2)"),
    Line2D([0], [0], marker="s", linestyle="",
           markerfacecolor="#bdb978", markeredgecolor="#8a8055",
           markersize=10, label=r"Background (non-moss surface, UAV NDVI $\leq$ 0.2)"),
]


# ---------------------------------------------------------------------------
# D -- Total moss area per UAV ortho (m^2)
# ---------------------------------------------------------------------------
MCD43A4_PIXEL_AREA = 500.0 * 500.0  # m^2 (used in panel F below)
colors = [SITE_COLORS[s] for s in SITES]
xpos = np.arange(len(SITES))

areas = [summary[s]["moss_area_m2"] for s in SITES]

ax_D.bar(
    xpos, areas, color=colors, edgecolor="black",
    linewidth=0.7, width=0.62,
)
for x, a in zip(xpos, areas):
    ax_D.text(x, a + max(areas) * 0.02, f"{a:.0f}",
              ha="center", va="bottom", fontsize=9.5, color="0.15")

ax_D.set_xticks(xpos)
ax_D.set_xticklabels([SITE_DISPLAY[s] for s in SITES])
ax_D.set_ylabel(r"Moss area in UAV ortho (m$^2$)")
ax_D.set_title("D. Total moss area per UAV ortho")
ax_D.set_ylim(0, max(areas) * 1.20)
ax_D.grid(axis="y", color="0.85", linewidth=0.5, zorder=0)
ax_D.set_axisbelow(True)


# ---------------------------------------------------------------------------
# E -- sub-pixel fraction of MCD43A4 pixel
# ---------------------------------------------------------------------------
pct = [summary[s]["moss_pct_pixel"] for s in SITES]
ax_E.bar(
    xpos, pct, color=colors, edgecolor="black", linewidth=0.7, width=0.62,
)
for x, p in zip(xpos, pct):
    ax_E.text(x, p + 0.002, f"{p:.4f}%",
              ha="center", va="bottom", fontsize=9.5, color="0.15")
ax_E.set_xticks(xpos)
ax_E.set_xticklabels([SITE_DISPLAY[s] for s in SITES])
ax_E.set_ylabel("Moss fraction of one MCD43A4 pixel (%)")
ax_E.set_title("E. Sub-pixel moss fraction")
ax_E.set_ylim(0, max(pct) * 1.22)
ax_E.grid(axis="y", color="0.85", linewidth=0.5, zorder=0)
ax_E.set_axisbelow(True)


# ---------------------------------------------------------------------------
# F -- UAV footprint vs MCD43A4 500-m pixel, drawn to scale
# ---------------------------------------------------------------------------
PIX = 500.0
# Outline of the MCD43A4 pixel; explicit label registered so it joins the
# legend below the panel.
mcd_patch = Rectangle(
    (0, 0), PIX, PIX,
    facecolor="white", edgecolor="black", linewidth=1.2,
    label="MCD43A4 pixel (500 m × 500 m)",
)
ax_F.add_patch(mcd_patch)
patch_handles = [mcd_patch]

# UAV rectangles laid out horizontally in the lower-left interior of the
# pixel. Each box is labeled inline with its area fraction of one MCD43A4
# pixel; site name and physical dimensions go in the legend below.
anchor_x, anchor_y = 30, 30
gap = 18
for site in SITES:
    w = summary[site]["uav_w"]
    h = summary[site]["uav_h"]
    frac = summary[site]["valid_area_m2"] / MCD43A4_PIXEL_AREA * 100.0
    p = Rectangle(
        (anchor_x, anchor_y), w, h,
        facecolor=SITE_COLORS[site], alpha=0.65,
        edgecolor="black", linewidth=0.6,
        label=f"{SITE_DISPLAY[site]} ({w:.0f} m × {h:.0f} m)",
    )
    ax_F.add_patch(p)
    patch_handles.append(p)
    # Inline area-fraction label, inside the box (centered)
    ax_F.text(
        anchor_x + w / 2.0, anchor_y + h / 2.0,
        f"{frac:.2f}%",
        ha="center", va="center",
        fontsize=9.5, color="black", fontweight="bold",
    )
    anchor_x += w + gap

ax_F.set_xlim(-20, PIX + 20)
ax_F.set_ylim(-20, PIX + 20)
ax_F.set_aspect("equal")
ax_F.set_xlabel("Distance (m)")
ax_F.set_ylabel("Distance (m)")
ax_F.set_title("F. UAV footprint vs MCD43A4 pixel scale")
ax_F.set_xticks([0, 100, 200, 300, 400, 500])
ax_F.set_yticks([0, 100, 200, 300, 400, 500])

# Compact 2x2 legend placed BELOW panel F, below the x-axis label, so it
# never overlaps the rectangles inside the panel.
ax_F.legend(
    handles=patch_handles,
    loc="upper center",
    bbox_to_anchor=(0.5, -0.20),
    ncol=2, frameon=False, fontsize=8.0,
    handlelength=1.3, handletextpad=0.5,
    columnspacing=1.4, labelspacing=0.4,
)


# Per-panel stats text in FIGURE coordinates, aligned at a single y just
# below the lowest map's bottom edge (panel C is the tallest map and
# therefore has the lowest bottom). This keeps the three text blocks on a
# common baseline so C's stats can't push down into panel F.
posA = ax_maps[0].get_position()
posC_bottom = min(ax.get_position().y0 for ax in ax_maps)
stats_y = posC_bottom - 0.010

for ax, site in zip(ax_maps, SITES):
    info = summary[site]
    pos = ax.get_position()
    cx = pos.x0 + pos.width / 2.0
    fig.text(
        cx, stats_y,
        f"Moss area = {info['moss_area_m2']:.0f} m$^2$\n"
        f"UAV cover = {info['moss_pct_uav']:.2f}%\n"
        f"MCD43A4 px cover = {info['moss_pct_pixel']:.3f}%",
        ha="center", va="top",
        fontsize=8.5, color="0.15", linespacing=1.35,
    )

# Shared legend rendered as a single boxed "special" panel. Placed in the
# whitespace immediately below panel A but ABOVE the common stats baseline,
# since A is the shortest map and that gap is otherwise unused. This keeps
# the legend clear of panel D below and of C's stats text.
legend_x = posA.x0 + posA.width / 2.0
legend_y = posA.y0 - 0.045  # top of legend, sitting in the gap below A
fig.legend(
    handles=legend_handles,
    loc="upper center",
    bbox_to_anchor=(legend_x, legend_y),
    ncol=1, frameon=True, fontsize=8.5,
    handlelength=1.2, handletextpad=0.5, labelspacing=0.45,
    borderpad=0.55,
    fancybox=True, edgecolor="0.6", facecolor="white",
)

fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
fig.savefig(OUT_PDF, bbox_inches="tight")
print(f"saved: {OUT_PNG}")
print(f"saved: {OUT_PDF}")
