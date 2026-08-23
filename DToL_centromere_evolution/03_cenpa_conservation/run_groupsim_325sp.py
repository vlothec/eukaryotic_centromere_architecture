#!/usr/bin/env python3
"""
Prepare inputs and run GroupSim on the bnni tree alignment,
comparing CENPA (420 seqs) vs H3-like (907 seqs).

Run from:  PhylogeneticProfiling/
  python3 19_curated_tree/run_groupsim_bnni.py
"""

import subprocess
import sys
from pathlib import Path
from collections import Counter

# ── paths ─────────────────────────────────────────────────────────────────────
BASE      = Path(__file__).parent
GROUPSIM  = Path("/home/jg2070/Desktop/groupsim-py3/src/groupsim.py")

ALN          = BASE / "cenpa430_H3_archaea10.aligned.clipkit.325sp.fasta"
H3_ALL_FASTA = BASE / "H3_all.aligned.fasta"

OUT_DIR = BASE / "split_entropy" / "groupsim"
OUT_DIR.mkdir(parents=True, exist_ok=True)

ARCHAEA_IDS = {
    "OLS22332.1","OLS24873.1","OLS21974.1","KKK41979.1","KXH71038.1",
    "OLS18261.1","OLS16336.1","BAD86478.1","OIO61677.1","OIO41945.1",
}

GAP = {"-"}

# ── helpers (same as split_entropy script) ────────────────────────────────────
def read_fasta(path):
    seqs, cur, buf = {}, None, []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line: continue
        if line.startswith(">"):
            if cur is not None: seqs[cur] = "".join(buf)
            cur = line[1:].split()[0]; buf = []
        else:
            buf.append(line)
    if cur is not None: seqs[cur] = "".join(buf)
    return seqs

def trim_columns(seqs, gap_threshold=0.85):
    ids  = list(seqs.keys())
    arr  = [seqs[i] for i in ids]
    n, L = len(arr), len(arr[0])
    keep = [i for i in range(L)
            if sum(1 for s in arr if s[i] in GAP) / n <= gap_threshold]
    trimmed = {sid: "".join(arr[k][c] for c in keep) for k, sid in enumerate(ids)}
    return trimmed, keep

def write_fasta(seqs, path):
    with open(path, "w") as fh:
        for sid, seq in seqs.items():
            fh.write(f">{sid}\n{seq}\n")

# ── main ──────────────────────────────────────────────────────────────────────
def run_for_gap(gap):
    tag = str(gap).replace(".", "")
    print(f"\n── Gap threshold {gap} ──")

    # 1. Read & filter
    seqs_all = read_fasta(ALN)
    seqs     = {k: v for k, v in seqs_all.items() if k not in ARCHAEA_IDS}

    # Read curated groups from groups_gap{tag}.txt (422 CENPA / 897 H3)
    curated_grp_file = OUT_DIR / f"groups_gap{tag}.txt"
    if curated_grp_file.exists():
        groups = {"CENPA": [], "H3": []}
        for line in curated_grp_file.read_text().splitlines():
            grp, ids = line.strip().split(":", 1)
            groups[grp] = [x for x in ids.split(",") if x in seqs]
    else:
        h3_pool = set(read_fasta(H3_ALL_FASTA).keys())
        groups = {"CENPA": [], "H3": []}
        for sid in seqs:
            groups["H3" if sid in h3_pool else "CENPA"].append(sid)
    print(f"  CENPA: {len(groups['CENPA'])}  H3: {len(groups['H3'])}")

    # 2. Trim
    trimmed, keep = trim_columns(seqs, gap_threshold=gap)
    print(f"  columns kept: {len(keep)}")

    # 3. Write trimmed FASTA (preserving group order: CENPA first, then H3)
    ordered = {sid: trimmed[sid] for sid in groups["CENPA"] + groups["H3"]
               if sid in trimmed}
    aln_out = OUT_DIR / f"trimmed_cenpa_h3_gap{tag}.fasta"
    write_fasta(ordered, aln_out)
    print(f"  alignment written: {aln_out.name}")

    # 4. Write group definition file
    grp_file = OUT_DIR / f"groups_gap{tag}.txt"
    with open(grp_file, "w") as fh:
        fh.write("CENPA:" + ",".join(groups["CENPA"]) + "\n")
        fh.write("H3:"    + ",".join(groups["H3"])    + "\n")
    print(f"  group file written: {grp_file.name}")

    # 5. Run GroupSim
    out_prefix = str(OUT_DIR / f"groupsim_gap{tag}")
    cmd = [
        sys.executable, str(GROUPSIM),
        "-k", str(grp_file),
        "-o", out_prefix,
        "-c", str(gap),       # column gap cutoff = same threshold as trim
        "-g", "0.5",          # group gap cutoff (relax slightly — 50% within group)
        "-w", "3",            # conservation window
        "-l", "0.7",          # lambda
        "-n",                 # map scores to [0,1]
        str(aln_out)
    ]
    print(f"  running GroupSim …")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("STDERR:", result.stderr[-2000:])
        raise RuntimeError(f"GroupSim failed (gap={gap})")
    print(result.stderr.strip())
    print(f"  done → {out_prefix}.*")

if __name__ == "__main__":
    for gap in [0.80, 0.85, 0.90]:
        run_for_gap(gap)
    print("\nAll done. Outputs in:", OUT_DIR)
