# Final WGCNA Repository Audit

## Outcome

The public-release staging repository now contains only the final revised GSE75665 WGCNA based on beta = 30 and 4,289 retained genes. No analysis was rerun, no result was recalculated, and no manuscript file was modified.

## Included files

### Scripts

- `analysis/WGCNA/scripts/03_GSE75665_signed_WGCNA.R`
- `analysis/WGCNA/scripts/01_build_revised_WGCNA_figure_and_tables.R`

Both files are exact copies of the archived final scripts.

### Inputs

- `analysis/WGCNA/data/processed/GSE75665_WGCNA_input_FPKM.csv`
- `analysis/WGCNA/data/metadata/GSE75665_subject_pairing.csv`
- `analysis/WGCNA/data/reference/Codonopsis_54_overlap_targets.csv`

The FPKM matrix was included because it is a compact processed matrix derived from the public GEO dataset GSE75665 and is necessary to reproduce the exact retained-gene universe. The repository documents the unrecovered upstream GeneID-to-symbol conversion rather than inferring or recreating it.

### Required final results

- `WGCNA_analysis_config.csv`
- `WGCNA_gene_filter_summary.csv`
- `WGCNA_sample_QC.csv`
- `WGCNA_soft_threshold_candidates.csv`
- `WGCNA_soft_threshold_decision.csv`
- `WGCNA_module_assignments.csv`
- `WGCNA_module_sizes.csv`
- `WGCNA_module_trait_mixed_effects.csv`
- `WGCNA_gene_module_metrics.csv`
- `WGCNA_candidate_gene_metrics.csv`
- `WGCNA_Codonopsis_target_enrichment.csv`

### Additional final provenance/results

- `WGCNA_module_eigengenes.csv`
- `WGCNA_primary_network.rds`
- `WGCNA_hub_genes.csv`
- `WGCNA_sample_outlier_decision.txt`
- `WGCNA_sessionInfo.txt`

## Final-version assertions checked

- retained genes = 4,289;
- MAD threshold = 0.50484;
- network type = signed;
- TOM type = signed;
- selected beta = 30;
- R2 threshold 0.85 was not reached;
- beta 30 selected by the documented negative-slope/nondegenerate-connectivity fallback rule;
- `minModuleSize = 20`;
- `deepSplit = 2`;
- `mergeCutHeight = 0.25`;
- three biological modules plus grey;
- subject-aware mixed model;
- global BH correction across nine module-effect tests;
- no module-effect association significant after BH correction.

## Legacy exclusion check

The repository contains no filename or result table identified as a legacy beta = 16 network, no 2,810-gene module assignment, no 12-module output, no legacy binary-trait module-correlation table, and no legacy “MMP9 hub” result. S4 CSV mirrors were omitted to avoid redundant copies with unclear precedence. The beta = 16 row remains in the complete soft-threshold candidate diagnostics because it is one prespecified candidate power among powers 1-30; it is not selected and is not a beta = 16 network result.

## Missing or constrained items

1. The final analysis script retains the archived project-relative path layout; an external user must map the three included inputs to those paths or change path variables only.
2. The final figure/table script contains an archived absolute project root and expects the master `Supplementary Tables.xlsx` workbook, which is outside this WGCNA-only module. The authoritative WGCNA result files needed for verification are included.
3. A fully automated script converting the upstream GEO GeneID FPKM table to the final gene-symbol FPKM matrix was not recovered. The exact final processed matrix is included with checksums.
4. No repository-wide `renv.lock` was available. Exact package versions are provided in `WGCNA_sessionInfo.txt`.
5. The persistent external host/DOI is not created by this staging operation.

These limitations concern path portability, upstream matrix provenance, and repository deployment; they do not create competing WGCNA result versions.
