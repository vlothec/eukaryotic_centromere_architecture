# 02 — Centromere-architecture ancestral-state reconstruction

Markov (Mk) models of centromere-architecture evolution across the calibrated 325-species
tree, testing the **α–ω cyclical hypothesis** (that centromere state can cycle between
satellite-based [α] and transposon-based [ω] architectures). Uses `phytools::fitMk` and
marginal ASR via `phytools::ancr`.

States: **H** (holocentric), **Sat** (satellite/α), **Trans** (transposon/ω),
**Mixed** (satellite+transposon on different chromosomes).

## Run order
```bash
bash run_pipeline.sh          # runs the steps below in order
```
| Script | Does |
|---|---|
| `00_prepare_inputs.R` | build per-dataset tree + tip-state annotation |
| `02_custom_models.R` | design matrices, incl. `ARD_irrevH` (H as an irreversible sink; direct Sat↔Trans allowed) |
| `37_test_cyclical_ARD_irrevH.R` | fit models; LRT (cyclical vs terminal-transposon); stochastic maps |
| `38_find_reversals.R` | detect X→Y→X state reversals on the mapped histories |
| `39_independent_cycles.R` | collapse reversals sharing an entry node into independent cycle events |
| `40_all_models_table.R` | AICc model-comparison table (ER, SYM, ARD, ARD_irrevH, TerminalTrans) |
| `43_cycles_tree_mk_parsimony.R` | ASR trees (marginal Mk + Fitch parsimony), cycles annotated |

## Key output
Winning model (`ARD_irrevH`) AICc table, ancestral-state trees, and independent-cycle
counts supporting the α–ω hypothesis.

## Inputs (not included)
- Calibrated tree from `01_species_tree_calibration/`.
- `branch_symbol_anno.tsv` (per-species architecture annotation) — in the full repository.
