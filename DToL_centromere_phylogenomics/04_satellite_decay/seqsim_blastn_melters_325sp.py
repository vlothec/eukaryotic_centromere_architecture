#!/usr/bin/env python3
"""
Satellite sequence similarity decay — Melters et al. (PMC4053949) approach,
implemented with BLASTN (local alignment) instead of edlib.

Strategy:
  - Build one BLAST database per species (sampled_repeats.fasta sequences).
  - Query each species against every other species.
  - Per species pair: mean % identity and fraction of queries with a
    significant hit (e-value < EVALUE_THRESH) among the top HSPs.
  - Node-average by MRCA divergence time.
  - Fit asymptotic exponential with background fixed at 25%.
  - Plot: (A) % identity decay, (B) detection rate vs divergence time.

BLAST parameters (Melters et al.): -task blastn -word_size 8
E-value: note that e-values are database-size-dependent; we normalise
         by reporting the fraction of significant hits rather than raw e-values.

Outputs:
  seqsim_blastn_325sp.tsv
  seqsim_blastn_summary.tsv
  seqsim_blastn_325sp.{pdf,png}
"""

import re, tempfile, subprocess
from pathlib import Path
from collections import defaultdict
import multiprocessing as mp
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

ROOT    = Path("/home/jg2070/Desktop/dtol_review_August")
PUB     = ROOT / "DToL_phylogenomics_publication_325genomes"
BASE    = ROOT / "2026_trees/annotation_centromeres/repeat_similarity"
FASTA   = BASE / "sampled_repeats_1000.fasta"
TREE_F  = PUB / "01_species_tree/outputs/full_325sp_calibrated_correlatedlambda01.nwk"  # winning model: correlated, lambda=0.1
TAX_F   = ROOT / "2026_trees/annotation_centromeres/centromere_code_to_species.tsv"
SAT_TAB = ROOT / "2026_trees/annotation_centromeres/organized/plots/centromere_length/genome_vs_centromere_length_satellite_species_table.tsv"
FIG_DIR = PUB / "05_satellite_similarity/figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

BLASTN     = "/home/jg2070/miniforge3/envs/of3_env/bin/blastn"
MAKEDB     = "/home/jg2070/miniforge3/envs/of3_env/bin/makeblastdb"
EVALUE_THRESH = 1e-5
BACKGROUND = 25.0   # random nucleotide identity
NCORES     = min(8, mp.cpu_count())

# ── load sequences per species ────────────────────────────────────────────────
print("Loading sequences …")
sp_seqs = defaultdict(list)
cur, buf = None, []
for line in FASTA.read_text().splitlines():
    if line.startswith(">"):
        if cur and buf: sp_seqs[cur].append("".join(buf))
        cur, buf = line[1:].split("__")[0].strip(), []
    else:
        buf.append(line.strip())
if cur and buf: sp_seqs[cur].append("".join(buf))
sp_seqs = {k: [s for s in v if len(s) >= 50] for k, v in sp_seqs.items()}
REMAP = {
    "morcroc":   "morcorc",
    "iyarmarma": "iyambarma",
    "iynyspin":  "iynysspin",
    "jahalocto": "jrhalocto",
    "jamemvect": "janemvect",
    "drroscani": "roscan_s",
    "dhquebobu": "dhquerobu",
    "dcastglyc": "drastglyc",
    "lpschltabe": "lpschtabe",
    "lpschltriq": "lpschtriq",
}
sp_seqs = {REMAP.get(k, k): v for k, v in sp_seqs.items()}
print(f"  {len(sp_seqs)} species, {sum(len(v) for v in sp_seqs.values()):,} seqs")
print(f"  Applied remapping for {len(REMAP)} old assembly codes")

# ── taxonomy & tree ───────────────────────────────────────────────────────────
from Bio import Phylo
tr = Phylo.read(str(TREE_F), "newick")
for c in tr.get_terminals(): c.name = c.name.replace(".fa","")
tip_map = {re.sub(r"[0-9.].*$","", t.name.lower()): t.name
           for t in tr.get_terminals()}

def get_mya(a, b):
    tA, tB = tip_map.get(a), tip_map.get(b)
    if not tA or not tB or tA == tB: return None
    try: return tr.distance(tA, tB) / 2
    except: return None

tax = pd.read_csv(TAX_F, sep="\t")
tax["code"] = tax["fasta_base"].str.lower().str.replace(r"[0-9.].*$","",regex=True)
vert  = {"Actinopterygii","Aves","Mammalia","Reptilia","Amphibia","Chondrichthyes"}
virid = {"Algae","Bryophyta","Dicots","Monocots","Gymnosperms"}

def get_group(code):
    rows = tax[tax["code"]==code]
    if rows.empty: return None
    t = rows["taxa1"].iloc[0]
    if t in vert:  return "Vertebrates"
    if t in virid: return "Viridiplantae"
    if t == "Fungi": return "Fungi"
    return "Invertebrates"

# Being present in the repeat annotation (sampled_repeats.fasta) *is* the satellite
# criterion; the only additional requirement is a tree tip for MRCA divergence times.
species = sorted([s for s in sp_seqs if s in tip_map and sp_seqs[s]])
print(f"  {len(species)} satellite species in tree")

# BLAST identities are sequence-only (tree-independent); cache them so that re-dating
# the analysis with a different chronogram only recomputes MRCA times, not the BLAST.
CACHE = BASE / "seqsim_blastn_raw_identities.tsv"
if not CACHE.exists():
    # ── write per-species fasta files into a temp dir ─────────────────────────
    tmpdir = Path(tempfile.mkdtemp(prefix="blast_sat_"))
    print(f"  Writing fastas to {tmpdir} …")
    sp_fa = {}
    for sp in species:
        fa = tmpdir / f"{sp}.fa"
        with open(fa,"w") as fh:
            for i, seq in enumerate(sp_seqs[sp]):
                fh.write(f">{sp}__{i}\n{seq}\n")
        sp_fa[sp] = fa

    # doubled-repeat fasta for the DB (Melters: repeats vs duplicated repeats, so
    # phase-shifted / boundary-offset monomers can align across the tandem junction)
    sp_fa_db = {}
    for sp in species:
        fa = tmpdir / f"{sp}.dup.fa"
        with open(fa,"w") as fh:
            for i, seq in enumerate(sp_seqs[sp]):
                fh.write(f">{sp}__{i}\n{seq}{seq}\n")
        sp_fa_db[sp] = fa

    # make BLAST databases
    print(f"  Building {len(species)} BLAST databases …")
    def make_db(sp):
        subprocess.run([MAKEDB,"-in",str(sp_fa_db[sp]),"-dbtype","nucl",
                        "-out",str(tmpdir/sp),"-logfile","/dev/null"],
                       capture_output=True)
    with mp.Pool(NCORES) as pool:
        pool.map(make_db, species)

# ── per-pair BLASTN worker ────────────────────────────────────────────────────
def run_blast_pair(args):
    spA, spB = args
    out = subprocess.run(
        [BLASTN,
         "-query",  str(sp_fa[spA]),
         "-db",     str(tmpdir/spB),
         "-task",   "blastn",
         "-word_size","8",
         "-reward","1","-penalty","-1","-gapopen","2","-gapextend","2","-dust","no",
         "-outfmt", "6 qseqid sseqid pident length evalue bitscore qlen",
         "-evalue", "10",
         "-max_target_seqs","1",   # best hit per query
         "-num_threads","1"],
        capture_output=True, text=True
    )
    n_q = len(sp_seqs[spA])
    if not out.stdout.strip():
        # no hits at all — global identity = 25% (random background)
        return spA, spB, BACKGROUND, 0.0, 0

    seen = {}   # qid → (pident, aln_len, evalue, qlen)
    for line in out.stdout.splitlines():
        p = line.split("\t")
        if len(p) < 7: continue
        qid = p[0]
        if qid not in seen:
            seen[qid] = (float(p[2]), int(p[3]), float(p[4]), int(p[6]))

    # global % identity (Melters et al.):
    # aligned portion: pident × aln_len; unaligned: 25% × (qlen - aln_len)
    global_ids = []
    evals = []
    for qid, (pident, aln_len, evalue, qlen) in seen.items():
        aln_len = min(aln_len, qlen)
        glob = (pident * aln_len + BACKGROUND * (qlen - aln_len)) / qlen
        global_ids.append(glob)
        evals.append(evalue)

    # queries with no hit at all → global id = 25
    n_no_hit = n_q - len(seen)
    global_ids += [BACKGROUND] * n_no_hit

    frac_sig = sum(1 for e in evals if e < EVALUE_THRESH) / n_q
    return spA, spB, float(np.mean(global_ids)), frac_sig, len(seen)

import itertools
if CACHE.exists():
    print(f"\nReusing cached BLAST identities: {CACHE.name} (skipping the all-vs-all BLASTN)")
    _c = pd.read_csv(CACHE, sep="\t")
    results = list(_c[["spA","spB","mean_pct_id","frac_sig","n_hits"]].itertuples(index=False, name=None))
else:
    pairs = [(a, b) for a, b in itertools.combinations(species, 2)]
    print(f"\nRunning BLASTN: {len(pairs):,} pairs on {NCORES} cores …")
    with mp.Pool(NCORES) as pool:
        results = pool.map(run_blast_pair, pairs, chunksize=5)
    pd.DataFrame(results, columns=["spA","spB","mean_pct_id","frac_sig","n_hits"]).to_csv(CACHE, sep="\t", index=False)
    print(f"Cached raw BLAST identities -> {CACHE.name}")

# ── build dataframe ───────────────────────────────────────────────────────────
rows = []
for spA, spB, mean_id, frac_sig, n_hits in results:
    mya = get_mya(spA, spB)
    if mya is None or mya == 0: continue
    gA, gB = get_group(spA), get_group(spB)
    group = gA if gA == gB else "Cross-group"
    rows.append(dict(spA=spA, spB=spB, mya=mya, mean_pct_id=mean_id,
                     frac_sig=frac_sig, n_hits=n_hits, group=group))

df = pd.DataFrame(rows)
df.to_csv(BASE / "seqsim_blastn_melters_325sp.tsv", sep="\t", index=False)
print(f"Saved: seqsim_blastn_melters_325sp.tsv  ({len(df):,} pairs)")

# ── node averaging ─────────────────────────────────────────────────────────────
within = df[df["group"] != "Cross-group"].copy()
within["mya"] = within["mya"].round(3)   # collapse float path-sum noise: one MRCA node = one point
navg_id  = within.groupby(["mya","group"]).agg(
    sim=("mean_pct_id","mean"), n=("mean_pct_id","count")).reset_index()
navg_sig = within.groupby(["mya","group"]).agg(
    frac_sig=("frac_sig","mean"), n=("frac_sig","count")).reset_index()

# ── exponential fit (C free — empirical floor) ────────────────────────────────
def exp_free(t, A, lam, C): return A * np.exp(-lam * t) + C

pal = {"Vertebrates":"#1565C0","Invertebrates":"#E65100","Viridiplantae":"#2E7D32"}
groups = ["Vertebrates","Invertebrates","Viridiplantae"]
summary = []

fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharey="row")

for col, grp in enumerate(groups):
    col_p = pal[grp]

    # Panel A: % identity
    ax = axes[0][col]
    d  = navg_id[navg_id["group"]==grp].dropna(subset=["sim"])
    ax.scatter(d["mya"], d["sim"], s=d["n"]*0.8,
               color=col_p, alpha=0.55, zorder=3)
    hl_str = "n/a"
    try:
        d_fit = d[d["mya"] > 0]
        C0  = float(d_fit["sim"].quantile(0.10))
        A0  = float(d_fit["sim"].max()) - C0
        A_ub = A0 * 1.5
        popt,_ = curve_fit(exp_free, d_fit["mya"], d_fit["sim"],
                           p0=[A0, 0.02, C0],
                           bounds=([0, 1e-5, 0], [A_ub, 5, 100]),
                           maxfev=10000)
        A, lam, C = popt; hl = np.log(2)/lam; hl_str = f"{hl:.0f} My"
        t_pred = np.linspace(0, d["mya"].max(), 400)
        ax.plot(t_pred, exp_free(t_pred, A, lam, C), color=col_p, lw=1.8)
        ax.axhline(C, color=col_p, lw=0.7, ls=":", alpha=0.6)
        ax.text(0.97, 0.95, f"t½ = {hl_str}\nfloor = {C:.1f}%",
                transform=ax.transAxes,
                ha="right", va="top", fontsize=10, color=col_p, fontweight="bold")
        summary.append(dict(group=grp, A=A, lambda_=lam, C_empirical=C,
                            halflife_My=hl, n_node_avg=len(d_fit)))
        print(f"  {grp:15s}  t½={hl:.1f} My  A={A:.2f}  λ={lam:.5f}  C={C:.2f}")
    except Exception as e:
        print(f"  {grp}: fit failed — {e}")
    ax.set_title(grp, fontsize=11, fontweight="bold", color=col_p)
    ax.set_ylabel("Mean % identity" if col==0 else "", fontsize=10)
    ax.spines[["top","right"]].set_visible(False)
    ax.yaxis.grid(True, color="#f0f0f0")

    # Panel B: detection rate (fraction with significant BLAST hit)
    ax2 = axes[1][col]
    d2  = navg_sig[navg_sig["group"]==grp].dropna(subset=["frac_sig"])
    ax2.scatter(d2["mya"], d2["frac_sig"]*100, s=d2["n"]*0.8,
                color=col_p, alpha=0.55, zorder=3)
    ax2.axhline(5, color="#999", lw=0.8, ls="--", alpha=0.7)
    ax2.set_ylabel("Sequences with BLAST hit (%)\n(e-value < 1e-5)" if col==0 else "",
                   fontsize=9)
    ax2.set_xlabel("Divergence time (My)", fontsize=10)
    ax2.set_ylim(-2, 105)
    ax2.spines[["top","right"]].set_visible(False)
    ax2.yaxis.grid(True, color="#f0f0f0")

axes[0][0].set_ylim(0, 100)
fig.suptitle("Satellite sequence similarity decay — BLASTN all-vs-all, node-averaged (325-sp)\n"
             "Top: mean % identity (C free — empirical floor).  "
             "Bottom: detection rate (fraction with BLASTN hit, e < 1×10⁻⁵).",
             fontsize=11, fontweight="bold")
plt.tight_layout()
for ext in ("pdf","png"):
    fig.savefig(FIG_DIR/f"seqsim_blastn_melters_325sp.{ext}",
                dpi=300 if ext=="png" else None,
                bbox_inches="tight", facecolor="white")
    print(f"Saved: seqsim_blastn_325sp.{ext}")
plt.close()

pd.DataFrame(summary).to_csv(BASE/"seqsim_blastn_melters_summary.tsv", sep="\t", index=False)

if "tmpdir" in globals():
    import shutil; shutil.rmtree(tmpdir, ignore_errors=True)
print("Done.")
