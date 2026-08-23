# 01 — Species tree & time calibration (methods)

Verified against the FastSpeciesTree run
`…/2026_trees/fast_species_tree_iqtree_only325/Results/` and the calibration tables
in `data/` (`calibration_table_325sp_supp_64.tsv`, `over_calib.tsv`). Superscript
numbers are manuscript reference markers.

## Species-tree inference
A maximum-likelihood species tree was inferred from the 325 DToL proteomes using
FastSpeciesTree⁶⁰ (v1.0) in a single "sensitive" run, which performs ortholog
identification, alignment and tree inference end-to-end. Each proteome was searched with
DIAMOND¹¹¹ blastp against the eukaryota_odb10 BUSCO reference orthologs; for each
single-copy ortholog the best-hit local alignment was anchored to the reference-query
coordinates and gap-padded to the query length, and these per-ortholog alignments were
concatenated into a query-anchored pseudo-alignment (orthologs absent from a proteome
encoded as gaps). Orthologs present in at least 80% of species (no more than 20% missing
taxa) were retained, yielding **194 single-copy loci**, and species with excessive gaps
were excluded. After FastSpeciesTree's internal column trimming, the supermatrix
comprised **325 taxa × 52,754 amino-acid sites (4.22% missing data)**. Within the same
run, a partitioned maximum-likelihood tree was inferred with IQ-TREE¹¹³ (**v2.3.4**)
under an **edge-linked proportional partition model** (`-p`), with ModelFinder selecting
the best-fit substitution model independently per locus from LG, JTT, Q.INSECT, Q.YEAST,
Q.BIRD, Q.MAMMAL and Q.PLANT, and among +I, +G and +I+G rate models. Across the 194
partitions, Q.INSECT (171 loci) and LG (15) were selected most often, with Q.YEAST (6),
Q.PLANT (1) and JTT (1) elsewhere (Q.BIRD and Q.MAMMAL were offered but never selected).
Branch support was assessed with 1,000 ultrafast bootstrap replicates (UFBoot) and 1,000
SH-aLRT replicates.

FastSpeciesTree was invoked as (sensitive mode → IQ-TREE):
```bash
python FastSpeciesTree.py -f all_proteomes_only325/ -o fast_species_tree_iqtree_only325 \
  -s sensitive -t 32
```
which internally ran the partitioned IQ-TREE command:
```bash
iqtree -T 32 -s trim_psuedo_alignment.fasta -p trim_IQTree_Partition_file.partitions \
  -B 1000 --alrt 1000 -st AA \
  -mset LG,JTT,Q.BIRD,Q.MAMMAL,Q.INSECT,Q.PLANT,Q.YEAST -mrate I,G,I+G -m MFP
```

> **Note (no MAFFT/trimAl).** FastSpeciesTree builds a DIAMOND-anchored pseudo-alignment
> and trims columns itself; MAFFT and trimAl are *not* part of this pipeline (its env is
> diamond + iqtree + veryfasttree only). The env pins iqtree 3.0.1, but the actual run
> used the `iqtree` on PATH — **v2.3.4** per the run log.

## Time calibration
The tree was **midpoint-rooted**, and calibration nodes were placed at well-supported
clades, each defined by the most recent common ancestor of two sampled taxa and assigned
a minimum–maximum age constraint. **Sixty-two constraints** were used in total, roughly
proportional to taxon sampling: **one root constraint** (Eukaryota, 1,085–1,671 Mya) and
**61 internal nodes** distributed across **Opisthokonta (n=1), Metazoa (n=31),
Viridiplantae (n=21) and Fungi (n=8)**. The bounds were taken from the corresponding
TimeTree 5⁶³ divergence-time confidence intervals — retrieved automatically with PAReTT
(https://github.com/LSLeClercq/PAReTT) — for **58 nodes**; the remaining **4 nodes**,
which lacked a reported TimeTree range, were bounded by the TimeTree median ± 20%
(the *Patella*, *Phorcus*/*Steromphala*, *Geum* and *Leistus*/*Pterostichus* divergences).

The tree topology was time-calibrated using the `chronos` function of the **ape** R
package⁶². To choose a dating model, we compared a range of penalised-likelihood models
and benchmarked each against TimeTree by node-age concordance (**213 node-age comparisons
across 210 shared taxa**): correlated-, relaxed- and discrete-rates models each across the
smoothing parameter λ = 0, 0.1, 1 and 10, together with a strict clock and a tree
constrained only at the root. The **correlated-rates model with λ = 0.1** gave the closest
agreement with TimeTree node ages (Pearson r = 0.98, R² = 0.97) and was used for plotting
and all further analyses.

```r
chronos(tree, lambda = 0.1, model = "correlated", calibration = calib_df)
```
