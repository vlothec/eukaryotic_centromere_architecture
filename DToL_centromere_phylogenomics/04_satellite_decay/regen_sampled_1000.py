#!/usr/bin/env python3
"""Regenerate the representative satellite set for the 165 curated species
(all.sats master; excludes the bat, Lycopus, and 8 holocentric species):
a random 10% of each species' repeats (>=50 bp), capped at 1000/species.
Two-pass reservoir sampling."""
import re, random, collections, math
random.seed(42)
SAT="/home/jg2070/Desktop/dtol_review_August/2026_trees/all.satellites.txt"
OUT="/home/jg2070/Desktop/dtol_review_August/2026_trees/annotation_centromeres/repeat_similarity/sampled_repeats_1000.fasta"
CAP=1000; FRAC=0.10
# species excluded in the master all.sats (bat + Lycopus + 8 holocentrics)
DROP={'mrhisin','dalyceuro','iglabmino','ihaelacum','ihicepurc','iilimauri',
      'ilpienapi','ioisceleg','iuloevari','iupsogibb'}
def sp_of(f): return re.sub(r'\.\d+$','',f[14].strip('"')).lower()
tot=collections.Counter()
for line in open(SAT):
    f=line.split()
    if len(f)<15: continue
    sp=sp_of(f)
    if sp in DROP: continue
    s=f[11].strip('"')
    if s and s!="NA" and len(s)>=50: tot[sp]+=1
target={sp:min(CAP,max(1,math.ceil(FRAC*n))) for sp,n in tot.items()}
res=collections.defaultdict(list); seen=collections.Counter()
for line in open(SAT):
    f=line.split()
    if len(f)<15: continue
    sp=sp_of(f)
    if sp in DROP: continue
    s=f[11].strip('"')
    if not s or s=="NA" or len(s)<50: continue
    k=target[sp]; j=seen[sp]; seen[sp]+=1
    if len(res[sp])<k: res[sp].append(s)
    else:
        r=random.randint(0,j)
        if r<k: res[sp][r]=s
with open(OUT,"w") as o:
    n=0
    for sp,seqs in res.items():
        for i,s in enumerate(seqs): o.write(f">{sp}__{i}\n{s}\n"); n+=1
print(f"species={len(res)}  total seqs={n}  at cap 1000: {sum(1 for sp in res if len(res[sp])==CAP)}")
