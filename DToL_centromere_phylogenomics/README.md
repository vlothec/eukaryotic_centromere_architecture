# DToL centromere phylogenomics

Core analyses for the Darwin Tree of Life centromere-evolution study (325 species): a
time-calibrated species tree and, on it, centromere-architecture evolution, CENP-A/CENH3
conservation, and satellite-DNA similarity decay.

This is a condensed, publication-oriented subset — the essential scripts, a few key figures,
and the trees/calibration tables needed to follow the analysis. Each sub-directory has its own
**README** giving the run order and, for every step, its input → script → output. The full
analysis repository (all figures and intermediate tables) is at
https://github.com/jacgonisa/DToL_phylogenomics

| module | analysis | key result |
|---|---|---|
| [`01_species_tree_calibration`](01_species_tree_calibration/) | FastSpeciesTree topology + `chronos` calibration (correlated rates, λ = 0.1; 62 TimeTree constraints) | calibrated 325-species chronogram |
| [`02_architecture_ASR`](02_architecture_ASR/) | Mk models of centromere-architecture evolution (7 models compared); α–ω cyclical hypothesis | `ARD_irrevH` best; independent cycles |
| [`03_cenpa_conservation`](03_cenpa_conservation/) | CENP-A/CENH3 (422) vs H3-like (897) entropy + GroupSim specificity-determining positions | SDPs in the loop-1/α2 CENP-A targeting domain |
| [`04_satellite_decay`](04_satellite_decay/) | satellite similarity decay vs divergence time + half-life | per-clade decay curves |

## External tools (referenced, not included)
- **FastSpeciesTree** — species-tree inference (DIAMOND-anchored BUSCO pseudo-alignment → IQ-TREE).
- **groupsim-py3** — GroupSim specificity-determining positions: https://github.com/jacgonisa/groupsim-py3
- **entropia-msa** — per-column alignment entropy: https://github.com/jacgonisa/entropia-msa

## Provided data
`01_species_tree_calibration/data/` holds the uncalibrated ML tree, the final calibrated
chronogram (`full_325sp_calibrated_correlatedlambda01.nwk`), and the 62 TimeTree calibration
points (PAReTT-retrieved). Larger inputs (assemblies, `all.satellites.txt`, the CENP-A/H3
alignment) live outside the repo; each module README states what they are. Scripts keep their
original file paths — edit the paths at the top of each to run.
