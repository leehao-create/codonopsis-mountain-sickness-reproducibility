# Analysis modules

| Module | Authoritative content | Entry point / validation |
|---|---|---|
| Network pharmacology | Frozen 54-target list, Figure 1 PPI metrics and enrichment display values | `network_pharmacology/scripts/draw_Fig1_final_candidate.py` |
| Differential expression | GSE103940 manuscript-locked analysis; GSE75665 interaction/GSEA; GSE260910 manuscript HAPE analysis with post hoc confounding audit | numbered scripts under `differential_expression/scripts/` |
| WGCNA | Manuscript-primary GSE75665 2,810-gene/beta-16/12-module workflow; reviewer checks isolated separately | `WGCNA/manuscript_primary/scripts/01_GSE75665_manuscript_primary_WGCNA.R` |
| Machine learning | Frozen GSE103940 nested subject-level LOSO, stability and paired permutation analysis | scripts `01`-`04`; `99_validate_outputs.R` validates frozen tables |
| Single cell | Descriptive Human Cell Landscape blood/lung localization | `single_cell/scripts/01_single_cell_candidate_localization.R` |
| Cell composition | Expression-based lineage-marker proxy audit | `cell_composition/scripts/01_cell_composition_marker_audit.R` |
| Molecular docking | Final NR3C2/MMP9 reviewer-response rerun using 2AA2/1GKC | activate Conda environment, then `molecular_docking/scripts/run_all.sh` |

The `results/` directories are immutable frozen outputs. Reproduction scripts use module-local `reproduced_results/` or `reproduced_figures/` directories by default. The ML GSE75665 material is explicitly supplementary/exploratory. GSE260910 was used in the manuscript HAPE analysis, while the released post hoc audit documents complete disease-batch confounding; it is not included in the released supervised machine-learning workflow.

WGCNA materials under `WGCNA/reviewer_checks/` are reviewer-requested sensitivity or diagnostic checks and are not primary manuscript findings. The 4,289-gene/beta-30 exploration is preserved only under `WGCNA/archive/nonmanuscript_beta30/`.
