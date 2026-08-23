#!/usr/bin/env python3
"""
Split-entropy plots for the bnni tree alignment (420 CenpA vs 907 H3, no archaea).

Produces per gap threshold:
  A. split_entropy_bnni_pretty_gap*.{png,pdf}
       Combined-trim approach (original): columns trimmed using all 1327 sequences.
  B. split_entropy_bnni_pergroup_gap*.{png,pdf}
       Per-group masking: same shared columns, but entropy set to NaN at positions
       where the group's own gap fraction exceeds the threshold.
  C. msa_matrix_gap*.{png,pdf}
       MSA amino-acid colour matrix — CENP-A block on top, H3 block below.

Run from:  PhylogeneticProfiling/
  python3 19_curated_tree/split_entropy_bnni_tree.py
"""

import math
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
import matplotlib.gridspec as gridspec

# ── paths ─────────────────────────────────────────────────────────────────────
BASE   = Path(__file__).parent
PP_DIR = Path("/home/jg2070/Desktop/dtol_review_August/PhylogeneticProfiling")
ENTDIR = PP_DIR / "20_cenpa_analysis" / "cenpa_plus_H3_entropy"

ALN          = BASE / "cenpa430_H3_archaea10.aligned.clipkit.325sp.fasta"
H3_ALL_FASTA = BASE / "H3_all.aligned.fasta"
STRIDE_H3    = ENTDIR / "AF-P59226-F1.stride"
STRIDE_CENPA = ENTDIR / "AF-Q8RVQ9-F1.stride"

OUT_DIR = BASE / "split_entropy"
OUT_DIR.mkdir(parents=True, exist_ok=True)

MAPPED_H3_ID    = "drParJuda1_SUPER_1_000979.1"
MAPPED_CENPA_ID = "ddAraThal4_chr_1_000021.1"

ARCHAEA_IDS = {
    "OLS22332.1", "OLS24873.1", "OLS21974.1", "KKK41979.1", "KXH71038.1",
    "OLS18261.1", "OLS16336.1", "BAD86478.1", "OIO61677.1", "OIO41945.1",
}

GAP = {"-"}

# ── amino acid colour map (Taylor-inspired) ───────────────────────────────────
AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY-")
AA_COLORS = {
    "A": "#80a0f0", "C": "#f08080", "D": "#c048c0", "E": "#c048c0",
    "F": "#80a0f0", "G": "#f09048", "H": "#15a4a4", "I": "#80a0f0",
    "K": "#f01505", "L": "#80a0f0", "M": "#80a0f0", "N": "#15c015",
    "P": "#ffff00", "Q": "#15c015", "R": "#f01505", "S": "#15c015",
    "T": "#15c015", "V": "#80a0f0", "W": "#09412d", "Y": "#15a4a4",
    "-": "#f0f0f0",   # gap = light grey
}
AA_IDX   = {aa: i for i, aa in enumerate(AA_ORDER)}
CMAP_COLORS = [AA_COLORS[aa] for aa in AA_ORDER]
AA_CMAP  = ListedColormap(CMAP_COLORS, name="aa_taylor")


# ── helpers ───────────────────────────────────────────────────────────────────
def read_fasta(path: Path):
    seqs, cur, buf = {}, None, []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if cur is not None:
                seqs[cur] = "".join(buf)
            cur = line[1:].split()[0]
            buf = []
        else:
            buf.append(line)
    if cur is not None:
        seqs[cur] = "".join(buf)
    return seqs


def trim_columns(seqs, gap_threshold=0.85):
    ids  = list(seqs.keys())
    arr  = [seqs[i] for i in ids]
    n, L = len(arr), len(arr[0])
    keep = [i for i in range(L)
            if sum(1 for s in arr if s[i] in GAP) / n <= gap_threshold]
    trimmed = {sid: "".join(arr[k][c] for c in keep) for k, sid in enumerate(ids)}
    return trimmed, keep


def group_gap_fractions(group_seqs, keep):
    """Gap fraction per kept column, computed only within group_seqs."""
    ids = list(group_seqs.keys())
    arr = [group_seqs[i] for i in ids]
    n   = len(arr)
    fracs = np.array([
        sum(1 for s in arr if s[i] in GAP) / n
        for i in range(len(arr[0]))
    ])
    return fracs   # length = len(keep)


def entropy_per_col(seqs_dict, nan_mask=None):
    """Shannon entropy per column, normalised by log2(20).
    nan_mask: boolean array, True → set entropy to NaN at that position."""
    seqs    = list(seqs_dict.values())
    L       = len(seqs[0])
    max_ent = math.log2(20)
    ent     = np.zeros(L)
    for i in range(L):
        col = [s[i] for s in seqs if s[i] not in GAP]
        if not col:
            continue
        counts = Counter(col)
        total  = sum(counts.values())
        h = -sum((v / total) * math.log2(v / total) for v in counts.values())
        ent[i] = h / max_ent
    if nan_mask is not None:
        ent = ent.astype(float)
        ent[nan_mask] = np.nan
    return ent


def parse_stride_helix_ranges(path: Path):
    ranges = []
    for line in path.read_text().splitlines():
        if line.startswith("LOC  AlphaHelix"):
            parts = line.split()
            try:
                ranges.append((int(parts[3]), int(parts[6])))
            except Exception:
                continue
    return ranges


def seqpos_to_alnpos(aln_seq: str):
    mapping, pos = {}, 0
    for i, aa in enumerate(aln_seq, start=1):
        if aa != "-":
            pos += 1
            mapping[pos] = i
    return mapping


def helix_ranges_in_aln(aln_seq: str, ranges):
    m = seqpos_to_alnpos(aln_seq)
    return [(m[s], m[e]) for s, e in ranges if s in m and e in m]


def project_to_trimmed(s, e, keep):
    inside = [k for k in keep if s <= k <= e]
    if inside:
        return (inside[0], inside[-1])
    nearest = lambda x: min(keep, key=lambda k: abs(k - x))
    a, b = nearest(s), nearest(e)
    return (min(a, b), max(a, b))


def add_helix_bands(ax, helix_cenpa, helix_h3):
    for s, e in helix_cenpa:
        ax.axvspan(s, e, ymin=0.93, ymax=0.99, color="#C62828", alpha=0.30, lw=0)
    for s, e in helix_h3:
        ax.axvspan(s, e, ymin=0.86, ymax=0.92, color="#1565C0", alpha=0.30, lw=0)


# ── plot A: combined-trim (original approach) ─────────────────────────────────
def plot_combined(df, helix_h3, helix_cenpa, out_prefix, gap, n_cenpa, n_h3):
    fig, ax = plt.subplots(figsize=(14, 4))

    ax.plot(df["pos"], df["H3"],    color="#1565C0", linewidth=0.8,
            alpha=0.85, label=f"H3-like (n={n_h3})")
    ax.plot(df["pos"], df["CENPA"], color="#C62828", linewidth=0.9,
            alpha=0.90, label=f"CENP-A (n={n_cenpa})")

    add_helix_bands(ax, helix_cenpa, helix_h3)

    thresh = np.nanpercentile(df["CENPA"], 90)
    top    = df[df["CENPA"] >= thresh]
    ax.scatter(top["pos"], top["CENPA"], s=10, color="#7B0000",
               alpha=0.7, zorder=3, label="Top 10% CENP-A positions")

    patch_cenpa = mpatches.Patch(color="#C62828", alpha=0.4, label="CENP-A α-helices")
    patch_h3    = mpatches.Patch(color="#1565C0", alpha=0.4, label="H3 α-helices")
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles + [patch_cenpa, patch_h3],
              labels  + ["CENP-A α-helices", "H3 α-helices"],
              fontsize=8, loc="upper right", framealpha=0.85)

    ax.set_xlabel(f"Alignment position (combined trim, gap ≤ {gap*100:.0f}%)", fontsize=10)
    ax.set_ylabel("Normalised Shannon entropy", fontsize=10)
    ax.set_title(
        f"Split entropy — bnni tree alignment  [combined trim]\n"
        f"CENP-A (n={n_cenpa}) vs H3-like (n={n_h3}), archaea excluded",
        fontsize=11, fontweight="bold"
    )
    ax.set_ylim(0, 1.08)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    plt.savefig(out_prefix.with_suffix(".png"), dpi=200)
    plt.savefig(out_prefix.with_suffix(".pdf"))
    plt.close()
    print(f"  saved {out_prefix.with_suffix('.png').name}")


# ── plot B: per-group masking ─────────────────────────────────────────────────
def plot_pergroup(df_pg, helix_h3, helix_cenpa, out_prefix, gap, n_cenpa, n_h3):
    fig, ax = plt.subplots(figsize=(14, 4))

    ax.plot(df_pg["pos"], df_pg["H3"],    color="#1565C0", linewidth=0.8,
            alpha=0.85, label=f"H3-like (n={n_h3})")
    ax.plot(df_pg["pos"], df_pg["CENPA"], color="#C62828", linewidth=0.9,
            alpha=0.90, label=f"CENP-A (n={n_cenpa})")

    add_helix_bands(ax, helix_cenpa, helix_h3)

    thresh = np.nanpercentile(df_pg["CENPA"].dropna(), 90)
    top    = df_pg[df_pg["CENPA"] >= thresh]
    ax.scatter(top["pos"], top["CENPA"], s=10, color="#7B0000",
               alpha=0.7, zorder=3, label="Top 10% CENP-A positions")

    # shade masked (NaN) regions
    cenpa_nan = df_pg["CENPA"].isna()
    h3_nan    = df_pg["H3"].isna()
    for i, (cn, hn) in enumerate(zip(cenpa_nan, h3_nan)):
        pos = df_pg["pos"].iloc[i]
        if cn:
            ax.axvspan(pos - 0.5, pos + 0.5, color="#C62828", alpha=0.08, lw=0)
        if hn:
            ax.axvspan(pos - 0.5, pos + 0.5, color="#1565C0", alpha=0.08, lw=0)

    patch_cenpa  = mpatches.Patch(color="#C62828", alpha=0.4, label="CENP-A α-helices")
    patch_h3     = mpatches.Patch(color="#1565C0", alpha=0.4, label="H3 α-helices")
    patch_masked = mpatches.Patch(color="grey",    alpha=0.25, label="Masked (group gap > threshold)")
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles + [patch_cenpa, patch_h3, patch_masked],
              labels  + ["CENP-A α-helices", "H3 α-helices", "Masked (group gap > thr)"],
              fontsize=8, loc="upper right", framealpha=0.85)

    ax.set_xlabel(f"Alignment position (per-group mask, gap ≤ {gap*100:.0f}%)", fontsize=10)
    ax.set_ylabel("Normalised Shannon entropy", fontsize=10)
    ax.set_title(
        f"Split entropy — bnni tree alignment  [per-group masking]\n"
        f"CENP-A (n={n_cenpa}) vs H3-like (n={n_h3}), archaea excluded",
        fontsize=11, fontweight="bold"
    )
    ax.set_ylim(0, 1.08)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    plt.savefig(out_prefix.with_suffix(".png"), dpi=200)
    plt.savefig(out_prefix.with_suffix(".pdf"))
    plt.close()
    print(f"  saved {out_prefix.with_suffix('.png').name}")


# ── plot C: MSA matrix ────────────────────────────────────────────────────────
def plot_msa_matrix(trimmed, groups, keep, helix_h3_trim, helix_cenpa_trim,
                    out_prefix, gap, n_cenpa, n_h3):
    cenpa_ids = sorted(groups["CENPA"])
    h3_ids    = sorted(groups["H3"])
    ncols     = len(keep)

    def seqs_to_matrix(ids):
        mat = np.zeros((len(ids), ncols), dtype=np.int8)
        for r, sid in enumerate(ids):
            seq = trimmed.get(sid, "-" * ncols)
            for c, aa in enumerate(seq):
                mat[r, c] = AA_IDX.get(aa.upper(), AA_IDX["-"])
        return mat

    print(f"  building MSA matrices ({n_cenpa} + {n_h3} × {ncols}) …")
    mat_cenpa = seqs_to_matrix(cenpa_ids)
    mat_h3    = seqs_to_matrix(h3_ids)

    # figure: two image panels stacked, shared x-axis
    h_cenpa = max(2.0, n_cenpa * 0.012)
    h_h3    = max(2.0, n_h3    * 0.007)
    fig     = plt.figure(figsize=(16, h_cenpa + h_h3 + 1.5))
    gs      = gridspec.GridSpec(2, 1, height_ratios=[n_cenpa, n_h3],
                                hspace=0.08, figure=fig)

    # shared x limits with a little padding
    xlim = (-0.5, ncols - 0.5)

    for ax, mat, title, helix_ref in [
        (fig.add_subplot(gs[0]), mat_cenpa,
         f"CENP-A  (n={n_cenpa})", helix_cenpa_trim),
        (fig.add_subplot(gs[1]), mat_h3,
         f"H3-like  (n={n_h3})",  helix_h3_trim),
    ]:
        ax.imshow(mat, aspect="auto", interpolation="nearest",
                  cmap=AA_CMAP, vmin=0, vmax=len(AA_ORDER) - 1,
                  origin="upper")

        # helix bars as coloured spans above the image
        for s, e in helix_ref:
            ax.axvspan(s - 1, e - 1, color="#444444", alpha=0.18, lw=0)

        ax.set_xlim(xlim)
        ax.set_ylabel(title, fontsize=9, labelpad=4)
        ax.tick_params(axis="y", left=False, labelleft=False)
        ax.tick_params(axis="x", labelsize=8)
        ax.spines[["top", "right", "left"]].set_visible(False)

    # shared x label on bottom panel
    fig.axes[-1].set_xlabel(
        f"Alignment position (combined trim, gap ≤ {gap*100:.0f}%)", fontsize=10)

    # amino-acid colour legend (compact, outside right)
    legend_patches = [
        mpatches.Patch(color=AA_COLORS[aa], label=aa)
        for aa in AA_ORDER if aa != "-"
    ]
    legend_patches.append(mpatches.Patch(color=AA_COLORS["-"], label="gap"))
    fig.legend(handles=legend_patches,
               title="Amino acid", title_fontsize=8,
               fontsize=7, loc="center right",
               bbox_to_anchor=(1.0, 0.5), ncol=1,
               framealpha=0.9, borderpad=0.5)

    fig.suptitle(
        f"MSA colour matrix — bnni tree alignment  (gap ≤ {gap*100:.0f}%)\n"
        f"CENP-A (n={n_cenpa}) on top · H3-like (n={n_h3}) below · grey bands = α-helices",
        fontsize=11, fontweight="bold", y=1.01
    )
    plt.savefig(out_prefix.with_suffix(".png"), dpi=150, bbox_inches="tight")
    plt.savefig(out_prefix.with_suffix(".pdf"),           bbox_inches="tight")
    plt.close()
    print(f"  saved {out_prefix.with_suffix('.png').name}")


# ── main ──────────────────────────────────────────────────────────────────────
def main():
    print("Reading alignment …")
    seqs_all = read_fasta(ALN)
    print(f"  total seqs: {len(seqs_all)}")

    seqs = {k: v for k, v in seqs_all.items() if k not in ARCHAEA_IDS}
    print(f"  after removing archaea: {len(seqs)}")

    print("Identifying H3 vs CENP-A from curated groups file …")
    groups_file = OUT_DIR / "groupsim" / "groups_gap08.txt"
    groups = {"H3": set(), "CENPA": set()}
    if groups_file.exists():
        for line in groups_file.read_text().splitlines():
            grp, ids = line.strip().split(":", 1)
            if grp in groups:
                groups[grp] = set(ids.split(",")) & set(seqs.keys())
    else:
        # fallback: infer from H3_ALL_FASTA
        h3_pool = set(read_fasta(H3_ALL_FASTA).keys())
        for sid in seqs:
            groups["H3" if sid in h3_pool else "CENPA"].add(sid)
    n_cenpa, n_h3 = len(groups["CENPA"]), len(groups["H3"])
    print(f"  CENP-A: {n_cenpa}  |  H3: {n_h3}")

    with open(OUT_DIR / "group_counts.txt", "w") as fh:
        for g, ids in groups.items():
            fh.write(f"{g}\t{len(ids)}\n")

    # structural annotation
    h3_ranges    = parse_stride_helix_ranges(STRIDE_H3)
    cenpa_ranges = parse_stride_helix_ranges(STRIDE_CENPA)
    aln_h3    = seqs[MAPPED_H3_ID]
    aln_cenpa = seqs[MAPPED_CENPA_ID]
    helix_h3_aln    = helix_ranges_in_aln(aln_h3,    h3_ranges)
    helix_cenpa_aln = helix_ranges_in_aln(aln_cenpa, cenpa_ranges)

    for gap in [0.80, 0.85, 0.90]:
        tag = str(gap).replace(".", "")
        print(f"\n── Gap threshold {gap} ──")

        trimmed, keep = trim_columns(seqs, gap_threshold=gap)
        print(f"  columns kept: {len(keep)}")

        keep_map = {orig: i + 1 for i, orig in enumerate(keep)}
        def project(aln_ranges):
            out = []
            for s, e in aln_ranges:
                proj = project_to_trimmed(s, e, keep)
                if proj and proj[0] in keep_map and proj[1] in keep_map:
                    out.append((keep_map[proj[0]], keep_map[proj[1]]))
            return out

        helix_h3_trim    = project(helix_h3_aln)
        helix_cenpa_trim = project(helix_cenpa_aln)

        # ── A: combined trim ──────────────────────────────────────────────────
        cenpa_sub = {k: trimmed[k] for k in groups["CENPA"] if k in trimmed}
        h3_sub    = {k: trimmed[k] for k in groups["H3"]    if k in trimmed}

        ent_cenpa = entropy_per_col(cenpa_sub)
        ent_h3    = entropy_per_col(h3_sub)
        ent_all   = entropy_per_col(trimmed)

        df = pd.DataFrame({
            "pos":   np.arange(1, len(ent_all) + 1),
            "all":   ent_all,
            "CENPA": ent_cenpa,
            "H3":    ent_h3,
        })
        df.to_csv(OUT_DIR / f"split_entropy_bnni_gap{tag}.tsv", sep="\t", index=False)

        plot_combined(df, helix_h3_trim, helix_cenpa_trim,
                      OUT_DIR / f"split_entropy_bnni_pretty_gap{tag}",
                      gap, n_cenpa, n_h3)

        # ── B: per-group masking ──────────────────────────────────────────────
        gf_cenpa = group_gap_fractions(cenpa_sub, keep)
        gf_h3    = group_gap_fractions(h3_sub,    keep)
        mask_cenpa = gf_cenpa > gap
        mask_h3    = gf_h3    > gap
        print(f"  masked positions — CENP-A: {mask_cenpa.sum()}  H3: {mask_h3.sum()}")

        ent_cenpa_pg = entropy_per_col(cenpa_sub, nan_mask=mask_cenpa)
        ent_h3_pg    = entropy_per_col(h3_sub,    nan_mask=mask_h3)

        df_pg = pd.DataFrame({
            "pos":   np.arange(1, len(ent_all) + 1),
            "CENPA": ent_cenpa_pg,
            "H3":    ent_h3_pg,
        })
        df_pg.to_csv(OUT_DIR / f"split_entropy_bnni_pergroup_gap{tag}.tsv",
                     sep="\t", index=False)

        plot_pergroup(df_pg, helix_h3_trim, helix_cenpa_trim,
                      OUT_DIR / f"split_entropy_bnni_pergroup_gap{tag}",
                      gap, n_cenpa, n_h3)

        # ── C: MSA matrix ─────────────────────────────────────────────────────
        plot_msa_matrix(trimmed, groups, keep,
                        helix_h3_trim, helix_cenpa_trim,
                        OUT_DIR / f"msa_matrix_gap{tag}",
                        gap, n_cenpa, n_h3)

    print("\nDone. Outputs in:", OUT_DIR)


if __name__ == "__main__":
    main()
