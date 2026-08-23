#!/usr/bin/env python3
"""Satellite-similarity decay / half-life figure using the chronos-correlated
calibrated tree. Reuses the saved BLASTN pair tsv (no rerun); only the MRCA
divergence times (mya) are recomputed from the chronos-correlated chronogram."""
import re
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from Bio import Phylo

ROOT = Path("/home/jg2070/Desktop/dtol_review_August")
PUB  = ROOT/"DToL_phylogenomics_publication_325genomes"
TSV  = ROOT/"2026_trees/annotation_centromeres/repeat_similarity/seqsim_blastn_melters_325sp.tsv"
TREE = PUB/"01_species_tree/outputs/full_325sp_calibrated_correlatedlambda01.nwk"  # winning model: correlated, lambda=0.1
OUT  = PUB/"05_satellite_similarity/figures"

# ── recompute mya from the chronos-correlated tree ────────────────────────────
tr = Phylo.read(str(TREE), "newick")
for c in tr.get_terminals(): c.name = c.name.replace(".fa","")
tip_map = {re.sub(r"[0-9.].*$","", t.name.lower()): t.name for t in tr.get_terminals()}
def get_mya(a, b):
    tA, tB = tip_map.get(a), tip_map.get(b)
    if not tA or not tB or tA == tB: return None
    try: return tr.distance(tA, tB) / 2
    except: return None

df = pd.read_csv(TSV, sep="\t")
df["mya"] = [get_mya(a, b) for a, b in zip(df["spA"], df["spB"])]
df = df.dropna(subset=["mya"])
df = df[df["mya"] > 0]

within = df[df["group"].isin(["Vertebrates","Invertebrates","Viridiplantae"])].copy()
within["group"] = within["group"].replace({"Vertebrates": "Chordata"})  # match paper (Fig 2/3) labelling
within["mya"] = within["mya"].round(3)   # collapse float path-sum noise: one MRCA node = one point
navg = within.groupby(["mya","group"]).agg(sim=("mean_pct_id","mean"),
                                           n=("mean_pct_id","count")).reset_index()

def exp_free(t, A, lam, C): return A*np.exp(-lam*t)+C
pal = {"Chordata":"#F72485","Invertebrates":"#3F37C9","Viridiplantae":"#8AC827"}  # paper palette (Fig 2/3)
groups = ["Chordata","Invertebrates","Viridiplantae"]
THR = 60   # identity threshold (%) for the "oldest node still above" marker

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
summary = []
for col, grp in enumerate(groups):
    cp = pal[grp]; ax = axes[col]
    d = navg[navg["group"]==grp].dropna(subset=["sim"])
    ax.scatter(d["mya"], d["sim"], s=d["n"]*0.8, color=cp, alpha=0.55, zorder=3)
    hl = np.nan; C = np.nan
    try:
        C0 = float(d["sim"].quantile(0.10)); A0 = float(d["sim"].max())-C0
        popt,_ = curve_fit(exp_free, d["mya"], d["sim"], p0=[A0,0.02,C0],
                           bounds=([0,1e-5,0],[A0*1.5,5,100]), maxfev=10000)
        A,lam,C = popt; hl = np.log(2)/lam
        t = np.linspace(0, d["mya"].max(), 400)
        ax.plot(t, exp_free(t,A,lam,C), color=cp, lw=1.8)
        ax.axhline(C, color=cp, lw=0.7, ls=":", alpha=0.6)
        ax.text(0.97,0.95,f"t½ = {hl:.0f} My\nfloor = {C:.1f}%", transform=ax.transAxes,
                ha="right", va="top", fontsize=10, color=cp, fontweight="bold")
    except Exception as e:
        print(grp, "fit failed", e)
    # (no >=60% marker: driven by over-dated shallow nodes, e.g. Schoenoplectus
    #  dated 17.6 My in-tree but ~3.4 My in TimeTree)
    ax.set_title(grp, fontsize=11, fontweight="bold", color=cp)
    ax.set_ylabel("Mean % identity" if col==0 else "", fontsize=10)
    ax.set_xlabel("Divergence time (My)", fontsize=10)
    ax.spines[["top","right"]].set_visible(False); ax.yaxis.grid(True, color="#f0f0f0")
    summary.append(dict(group=grp, half_life_My=round(hl,1), floor_pct=round(C,1),
                        n_nodes=len(d)))

axes[0].set_ylim(0, 100)
fig.suptitle("Satellite similarity decay — BLASTN all-vs-all, node-averaged (325-sp)\n"
             "Divergence times: chronos correlated.",
             fontsize=11, fontweight="bold")
plt.tight_layout()
for ext in ("png","pdf"):
    fig.savefig(OUT/f"seqsim_halflife_chronos_correlated_325sp.{ext}",
                dpi=300 if ext=="png" else None, bbox_inches="tight", facecolor="white")
    print("Saved:", OUT/f"seqsim_halflife_chronos_correlated_325sp.{ext}")

sdf = pd.DataFrame(summary)
sdf.to_csv(OUT/"halflife_chronos_correlated_325sp.tsv", sep="\t", index=False)
print(sdf.to_string(index=False))
