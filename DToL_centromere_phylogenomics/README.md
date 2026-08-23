# DToL centromere phylogenomics (325 species)

Core analyses from the Darwin Tree of Life centromere-evolution study: a
time-calibrated 325-species eukaryotic phylogeny and, on it, the evolution of
centromere architecture, the centromeric histone CENP-A/CENH3, and centromeric
satellite DNA.

This is a **condensed, publication-only** copy of the four core analyses. The full
analysis repository — every figure, model variant, and intermediate table — lives at
**https://github.com/jacgonisa/DToL_phylogenomics**.

## Contents

| Sub-directory | Analysis | Key output |
|---|---|---|
| [`01_species_tree_calibration`](01_species_tree_calibration/) | FastSpeciesTree topology + `chronos` time-calibration (correlated rates, λ=0.1; 62 TimeTree constraints) | calibrated 325-sp chronogram |
| [`02_architecture_ASR`](02_architecture_ASR/) | Mk models of centromere-architecture evolution; the α–ω cyclical hypothesis | ASR trees, model-comparison table, cycle counts |
| [`03_cenpa_conservation`](03_cenpa_conservation/) | CENP-A/CENH3 vs H3 per-column entropy + GroupSim specificity-determining positions (SDPs) | entropy profile, SDP set |
| [`04_satellite_decay`](04_satellite_decay/) | BLASTN all-vs-all satellite similarity decay vs divergence time + half-life | decay curves, per-clade half-lives |

## External tools (referenced, not vendored)

- **FastSpeciesTree** — species-tree inference from proteomes (DIAMOND-anchored BUSCO
  pseudo-alignment → partitioned IQ-TREE). See `01_species_tree_calibration/METHODS.md`.
- **groupsim-py3** — GroupSim specificity-determining positions (Python 3 port of
  Capra & Singh 2008): https://github.com/jacgonisa/groupsim-py3
- **entropia-msa** — per-column multiple-sequence-alignment entropy:
  https://github.com/jacgonisa/entropia-msa
- Also used: `ape::chronos` and `phytools` (R); PAReTT for TimeTree confidence intervals;
  NCBI BLASTN (v2.16.0+); PAML/MCMCTree (Bayesian dating cross-check).

## Notes on reproducing

- Scripts retain the **original absolute paths** (`/home/jg2070/…`) from the analysis
  environment; adjust the `BASE`/input paths at the top of each script to your layout.
- **Large inputs** (genome assemblies, `all.satellites.txt`, the CENP-A/H3 alignments)
  are not included here — each sub-directory README states what they are and how the
  in-repo tables were produced. The full repository above holds the derived data tables.
- Each analysis targets the 325 published species (the calibrated species-tree tips).
