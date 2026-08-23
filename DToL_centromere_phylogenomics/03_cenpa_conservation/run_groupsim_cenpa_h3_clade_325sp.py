#!/usr/bin/env python3
"""
GroupSim with clade-level weighting — CENP-A vs H3.

Clade weights: each sequence is weighted by 1 / n_clade, where n_clade is the
number of sequences from the same broad taxonomic group (Insects, Vertebrates,
Viridiplantae, Fungi, Other_Invertebrates) within the same comparison group
(CENP-A or H3). Weights are normalised within each group so the mean = 1.

Outputs: split_entropy/groupsim_weighted/
  groupsim_cenpa_h3_clade_gap08.tsv
  groupsim_cenpa_h3_clade_gap08.txt
"""

import math, re
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd
from scipy.stats import zscore as sp_zscore

BASE         = Path(__file__).parent
ALN          = BASE / "cenpa430_H3_archaea10.aligned.clipkit.325sp.fasta"
H3_ALL_FASTA = BASE / "H3_all.aligned.fasta"
TAX_F        = Path("/home/jg2070/Desktop/dtol_review_August/2026_trees/annotation_centromeres/centromere_code_to_species.tsv")
OUT_DIR      = BASE / "split_entropy" / "groupsim_weighted"
OUT_DIR.mkdir(parents=True, exist_ok=True)

GAP  = '-'
AA20 = list("ACDEFGHIKLMNPQRSTVWY")

ARCHAEA_IDS = {
    "OLS22332.1","OLS24873.1","OLS21974.1","KKK41979.1","KXH71038.1",
    "OLS18261.1","OLS16336.1","BAD86478.1","OIO61677.1","OIO41945.1",
}

INSECTS       = {"Coleoptera","Diptera","Hymenoptera","Lepidoptera","Hemiptera",
                 "Ephemeroptera","Odonata","Trichoptera","Neuroptera","Dermaptera",
                 "Plecoptera","Psocodea","Blattodea"}
VERTEBRATES   = {"Actinopterygii","Aves","Mammalia","Reptilia","Amphibia","Chondrichthyes"}
VIRIDIPLANTAE = {"Algae","Bryophyta","Dicots","Monocots","Gymnosperms"}

def taxa1_to_broad(t):
    if t in INSECTS:       return "Insects"
    if t in VERTEBRATES:   return "Vertebrates"
    if t in VIRIDIPLANTAE: return "Viridiplantae"
    if t == "Fungi":       return "Fungi"
    return "Other_Invertebrates"

def seq_to_code(sid):
    base = sid.split("_")[0].lower()
    return re.sub(r"[0-9.].*$", "", base)

def read_fasta(path):
    seqs, cur, buf = {}, None, []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line: continue
        if line.startswith(">"):
            if cur is not None: seqs[cur] = "".join(buf)
            cur = line[1:].split()[0]; buf = []
        else: buf.append(line.upper())
    if cur is not None: seqs[cur] = "".join(buf)
    return seqs

def trim_columns(seqs, thr):
    ids = list(seqs.keys()); arr = [seqs[i] for i in ids]
    n, L = len(arr), len(arr[0])
    keep = [c for c in range(L) if sum(1 for s in arr if s[c]==GAP)/n <= thr]
    return {sid: "".join(arr[k][c] for c in keep) for k,sid in enumerate(ids)}, keep

def clade_weights(seq_ids, group_idx_list, tax_df):
    code_to_taxa1 = dict(zip(
        tax_df["fasta_base"].str.lower().str.replace(r"[0-9.].*$","",regex=True),
        tax_df["taxa1"]))
    broad_arr = np.array([taxa1_to_broad(code_to_taxa1.get(seq_to_code(sid),"Unknown"))
                          for sid in seq_ids])
    w = np.ones(len(seq_ids))
    for idx in group_idx_list:
        clades = broad_arr[idx]
        counts = Counter(clades)
        for i, c in zip(idx, clades):
            w[i] = 1.0 / counts[c]
        grp_w = w[idx]; w[idx] = grp_w / grp_w.mean()
    return w, broad_arr

def weighted_freqs(col, indices, weights, group_gap_cutoff=0.5):
    g_col = [col[i] for i in indices]
    if sum(1 for c in g_col if c==GAP)/len(g_col) > group_gap_cutoff: return None
    pairs = [(aa, weights[i]) for i,aa in zip(indices,g_col) if aa!=GAP]
    if not pairs: return None
    wsum = sum(wi for _,wi in pairs)
    if wsum == 0: return None
    f = {}
    for aa,wi in pairs: f[aa] = f.get(aa,0.0) + wi/wsum
    return f

def groupsim_clade_score(col, group_idx_list, weights,
                          col_gap_cutoff=0.80, group_gap_cutoff=0.5):
    n_gap = sum(1 for c in col if c==GAP)
    if n_gap/len(col) > col_gap_cutoff: return None
    gfreqs = [weighted_freqs(col, idx, weights, group_gap_cutoff) for idx in group_idx_list]
    if any(gf is None for gf in gfreqs): return None
    within  = [sum(f**2 for f in gf.values()) for gf in gfreqs]
    between = []
    for g1 in range(len(gfreqs)):
        for g2 in range(g1+1, len(gfreqs)):
            all_aa = set(gfreqs[g1])|set(gfreqs[g2])
            between.append(sum(gfreqs[g1].get(a,0)*gfreqs[g2].get(a,0) for a in all_aa))
    return (1-n_gap/len(col)) * (np.mean(within) - (np.mean(between) if between else 0))

def js_div(col):
    bg = [0.074,0.052,0.045,0.054,0.025,0.034,0.054,0.074,0.026,0.068,
          0.099,0.058,0.025,0.047,0.039,0.057,0.051,0.013,0.032,0.073]
    pc = 0.001
    cnt = Counter(aa for aa in col if aa!=GAP); n = sum(cnt.values()) or 1
    fc1 = [cnt.get(aa,0)/n+pc for aa in AA20]; s1=sum(fc1); fc1=[x/s1 for x in fc1]
    fc2 = [x+pc for x in bg]; s2=sum(fc2); fc2=[x/s2 for x in fc2]
    r   = [0.5*fc1[i]+0.5*fc2[i] for i in range(20)]
    d   = sum(fc1[i]*math.log2(fc1[i]/r[i]) if fc1[i]>0 and r[i]>0 else 0 for i in range(20))
    d  += sum(fc2[i]*math.log2(fc2[i]/r[i]) if fc2[i]>0 and r[i]>0 else 0 for i in range(20))
    return (1-sum(1 for c in col if c==GAP)/len(col))*d/2

def norm_01(vals):
    v=[x for x in vals if x is not None]
    if not v: return vals
    lo,hi=min(v),max(v); r=hi-lo
    return [None if x is None else (1.0 if r==0 else (x-lo)/r) for x in vals]

def window_correct(raw, cols, window=3, lam=0.7):
    cons=[js_div(col) if s is not None else None for s,col in zip(raw,cols)]
    cons_n=norm_01(cons); out=[]
    for i,s in enumerate(raw):
        if s is None: out.append(None); continue
        lo,hi=max(i-window,0),min(i+window+1,len(raw))
        terms=[cons_n[j] for j in range(lo,hi) if j!=i and cons_n[j] is not None]
        wc=sum(terms)/len(terms) if terms else 0
        out.append((1-lam)*wc+lam*s)
    return out

def run(gap=0.80):
    tag = str(gap).replace('.','')
    print(f"\n── Gap threshold {gap} ──")

    seqs_all = read_fasta(ALN)
    seqs = {k:v for k,v in seqs_all.items() if k not in ARCHAEA_IDS}

    tag = str(gap).replace(".", "")
    curated = BASE / "split_entropy" / "groupsim" / f"groups_gap{tag}.txt"
    if curated.exists():
        _grp = {}
        for line in curated.read_text().splitlines():
            g, ids = line.strip().split(":", 1)
            _grp[g] = [x for x in ids.split(",") if x in seqs]
        cenpa_ids = _grp["CENPA"]
        h3_ids    = _grp["H3"]
    else:
        h3_pool   = set(read_fasta(H3_ALL_FASTA).keys())
        cenpa_ids = [sid for sid in seqs if sid not in h3_pool]
        h3_ids    = [sid for sid in seqs if sid in h3_pool]
    print(f"  CENPA: {len(cenpa_ids)}  H3: {len(h3_ids)}")

    trimmed, _ = trim_columns(seqs, gap)
    L = len(next(iter(trimmed.values())))
    print(f"  Columns after trim: {L}")

    ids_ord  = [s for s in cenpa_ids+h3_ids if s in trimmed]
    seqlist  = [trimmed[sid] for sid in ids_ord]
    cenpa_in = [s for s in cenpa_ids if s in trimmed]
    h3_in    = [s for s in h3_ids    if s in trimmed]
    idx_cenpa = list(range(len(cenpa_in)))
    idx_h3    = list(range(len(cenpa_in), len(ids_ord)))
    group_idx = [idx_cenpa, idx_h3]

    tax_df = pd.read_csv(TAX_F, sep='\t')
    print("  Computing clade weights …")
    w, broad_arr = clade_weights(ids_ord, group_idx, tax_df)

    for grp,name in [(idx_cenpa,"CENPA"),(idx_h3,"H3")]:
        cnt = Counter(broad_arr[grp])
        print(f"  {name} clades: {dict(cnt)}")
    print(f"  Weight range: {w.min():.4f} – {w.max():.4f}")

    cols = [[seq[c] for seq in seqlist] for c in range(L)]
    raw  = [groupsim_clade_score(cols[c], group_idx, w,
                                  col_gap_cutoff=gap, group_gap_cutoff=0.5)
            for c in range(L)]
    corrected = window_correct(raw, cols, window=3, lam=0.7)
    scores_n  = norm_01(corrected)

    valid = [(i,s) for i,s in enumerate(scores_n) if s is not None]
    vals  = np.array([s for _,s in valid])
    zs    = sp_zscore(vals) if len(vals)>1 else np.zeros(len(vals))
    zs_map = {i:z for (i,_),z in zip(valid,zs)}

    rows = []
    for c in range(L):
        s=scores_n[c]; z=zs_map.get(c)
        rows.append({"pos":c+1,"groupsim_clade":s,"z_clade":z,
                     "sig_clade":int(z is not None and z>=2.0)})
    df = pd.DataFrame(rows)
    out = OUT_DIR / f"groupsim_cenpa_h3_clade_gap{tag}.tsv"
    df.to_csv(out, sep='\t', index=False, float_format='%.6f')
    print(f"  Sig positions (z≥2): {df['sig_clade'].sum()}")
    print(f"  Saved: {out.name}")
    return df

if __name__ == "__main__":
    run(gap=0.80)
    run(gap=0.85)
