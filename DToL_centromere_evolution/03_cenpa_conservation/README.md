# 03 — CENP-A/CENH3 conservation & specificity

Two complementary views of how the centromeric histone CENP-A/CENH3 differs from
canonical histone H3 across the 325 species: per-column **conservation entropy** and
**GroupSim** specificity-determining positions (SDPs).

## Entropy (uses [entropia-msa](https://github.com/jacgonisa/entropia-msa))
- `split_entropy_325sp.py` — per-column Shannon entropy of the CENP-A/CENH3 vs H3
  alignment, split by group; gaps treated as a 21st character. Reads the curated
  group assignment (CENP-A vs H3) file.
- `plot_entropy_325sp.R` — positional-entropy profile with STRIDE-derived α-helix bands.
- `plot_panel_D_nomasking_325sp.py` — split-entropy panel (no-masking, filled version).

## GroupSim (uses [groupsim-py3](https://github.com/jacgonisa/groupsim-py3))
- `run_groupsim_325sp.py` — standard (unweighted) GroupSim: within- vs between-group
  residue similarity, window = 3, λ = 0.7, gaps as a 21st character, scores in [0,1];
  columns trimmed at >85 % gaps. Comparisons: CENP-A vs H3, and satellite- vs
  transposon-state CENP-A.
- `run_groupsim_cenpa_h3_clade_325sp.py` — clade-weighted variant (each sequence weighted
  by 1/n_clade over broad taxonomic groups; mean 1 per group) controlling for uneven
  sampling. **SDPs = clade-weighted z ≥ 2.**
- `plot_groupsim_cenpa_h3_gap085.R` — publication GroupSim figure (unweighted + weighted
  tracks) with the CENP-A helix track.

## Key output
Positional-entropy profile and the set of significant SDPs distinguishing CENP-A/CENH3
from H3 (concentrated in the histone-fold loop-1/α2 CENP-A targeting domain).

## Inputs (not included)
CENP-A/CENH3 + H3 protein alignment and the curated group-assignment file (in the full
repository); STRIDE `.stride` files for the helix bands.
