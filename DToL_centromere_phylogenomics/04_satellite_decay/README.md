# 04 — Satellite similarity decay & half-life

How centromeric satellite sequence similarity decays with species divergence time,
following the approach of Melters et al. (2013).

## Pipeline
1. `regen_sampled_1000.py` — from the satellite annotation (`all.satellites.txt`), draw a
   random 10 % sub-sample of each species' satellite monomers (≥ 50 bp) by reservoir
   sampling, capped at 1,000 sequences per species (seed 42) → `sampled_repeats_1000.fasta`
   (148,794 sequences across the satellite-bearing species).
2. `seqsim_blastn_melters_325sp.py` — all-vs-all pairwise BLASTN (v2.16.0+) between species'
   satellite sets. Each database sequence is a head-to-tail tandem **dimer** of the monomer
   (so phase-shifted monomers align across the repeat junction); best hit per query;
   `-word_size 8 -reward 1 -penalty -1 -gapopen 2 -gapextend 2 -dust no`. A global percent
   identity is recovered by extending each local hit over the full query length, assigning
   the 25 % identity expected for random nucleotide sequence to unaligned regions.
3. `halflife_chronos_correlated_325sp.py` — MRCA divergence times from the calibrated
   325-species tree; pairwise identities node-averaged by MRCA age; fit
   `H(t) = A·exp(−λt) + C` (background floor `C` a free parameter, `scipy.optimize.curve_fit`).
   Half-life `t½ = ln(2)/λ` (time to decay halfway from the initial value to the floor).

## Key output
Per-clade satellite similarity-decay curves and half-lives (Chordata / Invertebrates /
Viridiplantae).

## Inputs (not included)
- `all.satellites.txt` — satellite-monomer annotation (large; outside the repo).
- Calibrated tree from `01_species_tree_calibration/`.
