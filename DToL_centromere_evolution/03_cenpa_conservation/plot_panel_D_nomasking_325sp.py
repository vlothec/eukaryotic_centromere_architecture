#!/usr/bin/env python3
"""
Panel D — split entropy, no masking, filled-area aesthetic.

Differences from the original panel_D_split_entropy_325sp:
  - Uses combined-trim TSV (no per-group NaN masking)
  - Fills under each entropy curve with semi-transparent colour
    so CENP-A / H3 overlap regions are visible
  - No masked-region shading annotations

Output: 03_entropy/figures/panel_D_split_entropy_325sp_nomasking.{pdf,png}
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

BASE    = Path(__file__).parent
PUB_DIR = BASE.parent
SP_DIR  = PUB_DIR / "04_cenpa_phylogeny" / "split_entropy"
FIG_DIR = BASE / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ── data ──────────────────────────────────────────────────────────────────────
df = pd.read_csv(SP_DIR / "split_entropy_bnni_gap085.tsv", sep="\t")
helices = pd.read_csv(SP_DIR / "helix_positions_gap085.tsv", sep="\t")

pos    = df["pos"].values
cenpa  = df["CENPA"].values
h3     = df["H3"].values

from pathlib import Path as _P
_grp = (_P(__file__).parent.parent / "04_cenpa_phylogeny/split_entropy/groupsim/groups_gap085.txt").read_text().splitlines()
n_cenpa = len([x for l in _grp if l.startswith("CENPA:") for x in l.split(":",1)[1].split(",")])
n_h3    = len([x for l in _grp if l.startswith("H3:")    for x in l.split(":",1)[1].split(",")])

C_CENPA = "#C62828"
C_H3    = "#1565C0"

# ── figure ────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(13, 4.5))

# Filled areas — draw H3 first so CENP-A overlap shows on top
ax.fill_between(pos, 0, h3,
                color=C_H3, alpha=0.18, linewidth=0, zorder=1)
ax.fill_between(pos, 0, cenpa,
                color=C_CENPA, alpha=0.22, linewidth=0, zorder=2)

# Lines on top
ax.plot(pos, h3,    color=C_H3,    linewidth=0.75, alpha=0.90,
        label=f"H3-like  (n={n_h3})", zorder=4)
ax.plot(pos, cenpa, color=C_CENPA, linewidth=0.85, alpha=0.95,
        label=f"CENP-A  (n={n_cenpa})", zorder=5)

# Top-10% CENP-A dots
thresh   = np.nanpercentile(cenpa, 90)
top_mask = cenpa >= thresh
ax.scatter(pos[top_mask], cenpa[top_mask],
           s=9, color="#7B0000", alpha=0.75, zorder=6,
           label="Top 10% CENP-A positions")

# ── Helix annotation strips (narrow bands near top of y-axis) ─────────────────
for _, row in helices[helices["histone"] == "CENPA"].iterrows():
    ax.axvspan(row["start"], row["end"],
               ymin=0.935, ymax=0.990,
               color=C_CENPA, alpha=0.45, linewidth=0, zorder=3)

for _, row in helices[helices["histone"] == "H3"].iterrows():
    ax.axvspan(row["start"], row["end"],
               ymin=0.875, ymax=0.930,
               color=C_H3, alpha=0.45, linewidth=0, zorder=3)

# ── Axes & decoration ─────────────────────────────────────────────────────────
ax.set_xlim(pos.min(), pos.max())
ax.set_ylim(0, 1.01)
ax.set_xlabel("Alignment position (combined trim, gap ≤ 85%)", fontsize=11)
ax.set_ylabel("Normalised Shannon entropy", fontsize=11)
ax.set_yticks([0, 0.25, 0.50, 0.75, 1.00])
ax.spines[["top", "right"]].set_visible(False)
ax.yaxis.grid(True, color="#eeeeee", zorder=0)

# Legend
patch_cenpa_h = mpatches.Patch(color=C_CENPA, alpha=0.45, label="CENP-A α-helices")
patch_h3_h    = mpatches.Patch(color=C_H3,    alpha=0.45, label="H3 α-helices")
handles, labels = ax.get_legend_handles_labels()
ax.legend(
    handles + [patch_cenpa_h, patch_h3_h],
    labels  + ["CENP-A α-helices", "H3 α-helices"],
    fontsize=9, loc="upper right", framealpha=0.90,
    edgecolor="#cccccc", handlelength=1.4
)

ax.set_title(
    "CENP-A vs H3 split entropy — bnni tree (325-sp)",
    fontsize=12, fontweight="bold"
)
ax.text(0.5, -0.13,
        "No per-group gap masking applied.  "
        "Filled areas: CENP-A (red) / H3 (blue) — overlap regions show mixing.",
        ha="center", va="top", fontsize=8, color="#616161",
        transform=ax.transAxes)

plt.tight_layout(rect=[0, 0.04, 1, 1])

for ext in ("pdf", "png"):
    out = FIG_DIR / f"panel_D_split_entropy_325sp_nomasking.{ext}"
    fig.savefig(str(out),
                dpi=300 if ext == "png" else None,
                bbox_inches="tight", facecolor="white")
    print(f"Saved: {out}")

plt.close(fig)
