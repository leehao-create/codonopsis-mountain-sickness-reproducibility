# Cell-composition marker audit

This module evaluates aggregate expression-based lineage-marker proxies in the frozen bulk expression objects. It does not estimate cell proportions and is not a formal deconvolution. The final publication figure and all marker-level and score-level tables are included. GSE260910 values are descriptive only because disease and sequencing batch are completely confounded.

Marker scores are expression-based lineage proxies and are not quantitative estimates of cell proportions. The marker panel covers T cells, B cells, NK cells, neutrophils, monocytes, and erythroid/reticulocyte lineages. Submission-ready CSV mirrors `Table_S4A` through `Table_S4E` are provided under `results/supplementary_tables/`.

For GSE103940, the paired count/TMM/voom object and all-gene tables in `data/dependencies/` are retained solely to reproduce this reviewer-requested marker-score audit. They are not the primary differential-expression result reported in the manuscript. The manuscript-locked primary result (2,782 DEGs; 118 increased and 2,664 decreased at high altitude) and its FPKM-based reproduction materials are isolated under `analysis/differential_expression/results/GSE103940_manuscript/`.
