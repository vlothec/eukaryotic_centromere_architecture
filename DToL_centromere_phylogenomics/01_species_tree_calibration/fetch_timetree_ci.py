#!/usr/bin/env python3
"""Fetch real TimeTree median + CI (precomputed_age / precomputed_ci_low/high)
for every calibration node, via NCBI taxid -> timetree.org pairwise JSON API.
Cached to disk so reruns are free. Output: calibration_nodes_timetree.tsv."""
import csv, re, json, time, sys
import urllib.request, urllib.parse
from pathlib import Path

BASE = Path("/home/jg2070/Desktop/dtol_review_August")
PUB  = BASE/"DToL_phylogenomics_publication_325genomes/01_species_tree"
QC   = PUB/"outputs/calibration_qc"
TAX  = BASE/"2026_trees/annotation_centromeres/centromere_code_to_species.tsv"
CALIB= QC/"calib64_constraints_for_panelA.tsv"
CACHE= QC/"timetree_api_cache.json"
OUT  = QC/"calibration_nodes_timetree.tsv"

cache = json.loads(CACHE.read_text()) if CACHE.exists() else {"taxid": {}, "pair": {}}
UA = {"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) Chrome/120 Safari/537.36"}

def _get(url, tries=5):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            return json.loads(urllib.request.urlopen(req, timeout=25).read())
        except urllib.error.HTTPError as e:
            if e.code == 429:
                time.sleep(2*(i+1)); continue
            raise
        except Exception:
            time.sleep(1.5); continue
    return None

def taxid(sp):
    if sp in cache["taxid"]: return cache["taxid"][sp]
    q = urllib.parse.quote(f'"{sp}"[Scientific Name]')
    d = _get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={q}&retmode=json")
    time.sleep(0.4)  # NCBI: <3 req/s
    tid = None
    if d:
        idl = d.get("esearchresult", {}).get("idlist", [])
        tid = idl[0] if idl else None
    cache["taxid"][sp] = tid
    return tid

def pair(a, b):
    ta, tb = taxid(a), taxid(b)
    key = f"{ta}/{tb}"
    if ta is None or tb is None:
        return None, None, None, "no_taxid"
    if key in cache["pair"] and "adj" in cache["pair"][key]:
        r = cache["pair"][key]; return r["med"], r["lo"], r["hi"], r["adj"], r["status"]
    d = _get(f"https://timetree.org/api/pairwise/{ta}/{tb}/json")
    time.sleep(0.3)
    med = lo = hi = adj = None; status = "no_data"
    if d:
        s = d.get("studies") if isinstance(d.get("studies"), dict) else {}
        # TimeTree web: "Median Time"=precomputed_age, "Range"=ci_low/high, "Adjusted Time"=adjusted_age
        med = s.get("precomputed_age") or d.get("sum_simple_mol_time")
        lo  = s.get("precomputed_ci_low"); hi = s.get("precomputed_ci_high")
        adj = s.get("adjusted_age")
        def f2(x):
            try: return float(x) if x is not None else None
            except: return None
        med, lo, hi, adj = f2(med), f2(lo), f2(hi), f2(adj)
        if adj == 0: adj = None    # TimeTree returns 0 when no adjusted time exists
        status = "ok" if med is not None else "no_data"
    cache["pair"][key] = {"med": med, "lo": lo, "hi": hi, "adj": adj, "status": status}
    return med, lo, hi, adj, status

# species lookup
sp = {}
for r in csv.DictReader(open(TAX), delimiter="\t"):
    sp[re.sub(r"[0-9.].*$","", r["fasta"].lower())] = f'{r["genus"]} {r["species"]}'
# NCBI taxonomy synonyms (name in our data -> name NCBI/TimeTree recognise)
ALIAS = {"Potentilla anserina": "Argentina anserina"}

rows = []
for r in csv.DictReader(open(CALIB), delimiter="\t"):
    a = sp[re.sub(r"[0-9.].*$","", r["tip_a"].lower())]
    b = sp[re.sub(r"[0-9.].*$","", r["tip_b"].lower())]
    med, lo, hi, adj, st = pair(ALIAS.get(a, a), ALIAS.get(b, b))
    source = "TimeTree pairwise API" if med is not None else "unresolved"
    rows.append(dict(label=r["label"], clade=r["clade"], taxonA=a, taxonB=b,
                     our_age=r["actual_age"], tt_median=med, tt_ci_low=lo,
                     tt_ci_high=hi, tt_adjusted=adj, status=st, source=source))
    print(f"{r['label']:<16} {a} | {b}  -> med={med} adj={adj} CI=[{lo},{hi}] ({st})", flush=True)
    CACHE.write_text(json.dumps(cache))  # checkpoint each node

with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)
n_ci = sum(1 for x in rows if x["tt_ci_low"] not in (None, "None"))
print(f"\nWrote {OUT}: {len(rows)} nodes, {n_ci} with real TimeTree CI")
