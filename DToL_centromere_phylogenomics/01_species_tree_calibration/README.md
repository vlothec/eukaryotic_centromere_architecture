# 01 — Species tree & time calibration

A 325-species eukaryotic tree, time-calibrated with secondary (TimeTree) constraints.
Full methods: [`METHODS.md`](METHODS.md).

## Pipeline
1. **Species tree** — inferred with **FastSpeciesTree** (external): a DIAMOND-anchored
   pseudo-alignment of `eukaryota_odb10` BUSCO orthologs → partitioned IQ-TREE. See
   `METHODS.md` for the exact command and model settings. (Topology only; not re-run here.)
2. `fetch_timetree_ci.py` — retrieve TimeTree 5 divergence-time confidence intervals for
   the calibration nodes (via PAReTT).
3. `make_calibration_table_64_325sp.R` — assemble the calibration constraint table
   (62 constraints: 1 root + 61 internal; 58 with a TimeTree CI, 4 as median ± 20%).
4. `calibrate_chronos_correlated_325sp.R` — time-calibrate the midpoint-rooted tree with
   `ape::chronos` (correlated-rates model, smoothing λ = 0.1 — the selected model).
5. `plot_calibration_combined_benchmark_publication_325sp.R` — model-selection benchmark
   (correlated/relaxed/discrete × λ = 0, 0.1, 1, 10) scored against TimeTree by MRCA
   node-age concordance (Extended Data figure).

## Key output
`full_325sp_calibrated_correlatedlambda01.nwk` — the calibrated chronogram used by all
downstream analyses (root ≈ 1475 Ma).

## Inputs (not included)
- FastSpeciesTree proteome set + its output tree.
- `over_calib.tsv` / calibration tables (produced by steps 2–3; in the full repository).
