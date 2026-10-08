# Manuscript-primary WGCNA

This directory contains the WGCNA workflow and frozen outputs corresponding to the submitted manuscript.

The analysis used 20 GSE75665 samples, removed genes with missing values or zero total expression, retained 2,810 genes using `MAD > 1.5`, and used beta = 16. Soft-threshold diagnostics were evaluated using a signed-network setting, while module construction used the original `blockwiseModules` workflow with a signed TOM. The original call did not explicitly set `networkType`; it is preserved without altering its behavior. Other fixed parameters were `minModuleSize = 20`, `mergeCutHeight = 0.05`, `reassignThreshold = 0`, and `maxBlockSize = 2810`. The frozen workflow produced 12 modules, including grey/unclassified genes.

At beta = 16, the signed soft-threshold diagnostic reported a scale-free topology fit of R2 = 0.0419; the conventional R2 >= 0.85 criterion was not reached. Beta 16 was retained as the manuscript-used value and is not presented as satisfying that criterion.

The primary archived module-trait matrices use four group indicators (`high_AMS`, `high_non-AMS`, `low_AMS`, and `low_non-AMS`). Reviewer-requested binary-AMS correlations and multiplicity checks are isolated under `../reviewer_checks/`. No adjusted module association or module-size-aware candidate-target enrichment is interpreted as statistically significant, and MMP9 is not interpreted as a high-connectivity hub.

## Files

- `scripts/00_original_WGCNA_code_as_archived.R.txt`: byte-preserved archived code, including two console-output lines that make it non-parseable as an R script.
- `scripts/01_GSE75665_manuscript_primary_WGCNA.R`: repository-portable transcription that starts from the included processed matrix; core filtering, signed soft-threshold diagnostic, and original `blockwiseModules` call are unchanged.
- `data/GSE75665_original_4group_metadata.csv`: exact four-level group metadata used by the archived module-trait matrices.
- `results/WGCNA_gene_filter_summary.csv`: frozen filtering counts and threshold.
- `results/WGCNA_analysis_config.csv`: manuscript-used parameters and the explicit distinction between diagnostic and construction settings.
- `results/WGCNA_soft_threshold_selected_power.csv`: reported beta-16 diagnostic value.
- `results/WGCNA_module_assignments.csv`: 2,810 genes and their numeric/color module assignments.
- `results/WGCNA_module_sizes.csv`: sizes of all 12 modules.
- `results/WGCNA_module_eigengenes.csv`: archived module eigengenes.
- `results/WGCNA_module_trait_4group_correlations.csv` and `WGCNA_module_trait_4group_pvalues.csv`: original 12 x 4 matrices.

The processed input matrix and sample metadata are shared under `../data/`. The archived upstream GeneID-to-symbol conversion was not fully scripted; the exact processed matrix and its checksum are retained.
