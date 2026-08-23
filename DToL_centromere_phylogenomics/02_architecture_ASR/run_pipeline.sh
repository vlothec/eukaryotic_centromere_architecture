#!/usr/bin/env bash
# run_pipeline.sh — Centromere ASR pipeline (ARD_irrevH cyclical model)
# Run from the directory CONTAINING ASR_March2026_327species/:
#   cd /path/to/annotation_centromeres
#   bash ASR_March2026_327species/run_pipeline.sh

set -euo pipefail

SCRIPTS="ASR_March2026_327species/scripts/plotting"

echo "=== Step 1: Model comparison + stochastic mapping (slow) ==="
Rscript "$SCRIPTS/37_test_cyclical_ARD_irrevH.R"

echo "=== Step 2: Find reversal tips ==="
Rscript "$SCRIPTS/38_find_reversals.R"

echo "=== Step 3: Count independent cycle events ==="
Rscript "$SCRIPTS/39_independent_cycles.R"

echo "=== Step 4: All-models AICc table ==="
Rscript "$SCRIPTS/40_all_models_table.R"

echo "=== Step 5: Q matrix plots ==="
Rscript "$SCRIPTS/41_Qmatrix_plots.R"

echo "=== Step 6: Tree plots (Mk + Parsimony ASR) ==="
Rscript "$SCRIPTS/43_cycles_tree_mk_parsimony.R"

echo "=== Step 7: Formatted model table ==="
Rscript "$SCRIPTS/44_model_table.R"

echo ""
echo "Pipeline complete. Outputs in ASR_March2026_327species/outputs/"
