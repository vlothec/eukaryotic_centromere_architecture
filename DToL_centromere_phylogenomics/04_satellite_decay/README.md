# 04 — Satellite similarity decay & half-life

How centromeric satellite sequence similarity decays with species divergence time,
following Melters et al. (2013).

## Steps (order · input → script → output)

1. `regen_sampled_1000.py`
   in: satellite-monomer annotation (`all.satellites.txt`) ·
   out: `sampled_repeats_1000.fasta` — a random 10 % of each species' monomers (≥ 50 bp) by
   reservoir sampling, capped at 1,000 per species (seed 42; 148,794 sequences).

2. `seqsim_blastn_melters_325sp.py`
   in: `sampled_repeats_1000.fasta` · out: `seqsim_blastn_melters_325sp.tsv` (pairwise identities).
   All-vs-all BLASTN (v2.16.0+): each database sequence is a head-to-tail tandem **dimer** of the
   monomer (so phase-shifted monomers align across the repeat junction); best hit per query;
   `-word_size 8 -reward 1 -penalty -1 -gapopen 2 -gapextend 2 -dust no`. A global percent identity
   is recovered by extending each local hit over the full query length, assigning the 25 % identity
   expected for random sequence to unaligned regions.

3. `halflife_chronos_correlated_325sp.py`
   in: the pairwise tsv + calibrated tree (module 01) · out: decay curves + per-clade half-lives.
   MRCA divergence times from the calibrated tree; identities node-averaged by MRCA age; fit
   `H(t) = A·exp(−λt) + C` (background floor `C` a free parameter, `scipy.optimize.curve_fit`);
   half-life `t½ = ln(2)/λ`.

## Key output
Per-clade satellite similarity-decay curves and half-lives (Chordata / Invertebrates / Viridiplantae).

![decay](figures/seqsim_halflife_chronos_correlated_325sp.png)
