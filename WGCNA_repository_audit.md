# WGCNA Repository Version Audit

## Outcome

The repository has one clearly identified manuscript-primary WGCNA: the original GSE75665 2,810-gene, beta-16 workflow producing 12 modules. No biological analysis was rerun during repository assembly.

## Primary workflow trace

- Input: `analysis/WGCNA/data/processed/GSE75665_WGCNA_input_FPKM.csv`.
- Preserved code: `analysis/WGCNA/manuscript_primary/scripts/01_GSE75665_manuscript_primary_WGCNA.R`.
- Missing/zero-total filter: 20,751 genes remained.
- Variability filter: `MAD > 1.5`; 2,810 genes remained.
- Soft-threshold diagnostics: explicit signed-network setting; beta = 16; R2 = 0.0419; R2 >= 0.85 not reached.
- Module construction: original `blockwiseModules` call, with `TOMType = "signed"` and no explicit `networkType` argument.
- Other parameters: Pearson correlation, `minModuleSize = 20`, `mergeCutHeight = 0.05`, `reassignThreshold = 0`, and `maxBlockSize = 2810`.
- Frozen module assignment: 12 modules including grey.

The repository does not claim that an explicitly signed adjacency generated the 12-module result.

## Reviewer-check boundary

`analysis/WGCNA/reviewer_checks/` contains:

- binary AMS correlations and BH adjustment across 12 modules (minimum adjusted P = 0.0607);
- the original 12 x 4 group-indicator matrices and global BH adjustment across 48 values (minimum adjusted P = 0.387);
- module-size-aware candidate-target enrichment (all BH-adjusted P > 0.05);
- MMP9 diagnostics: turquoise, kME = 0.852, traditional GS = 0.4423, GS P = 0.0508, signed-adjacency kWithin = 226.2, rank 967/1,147.

These checks do not redefine the primary network. MMP9 is exploratory context and not a high-connectivity hub.

## Archived nonmanuscript analysis

The former 4,289-gene/beta-30 explicitly signed reanalysis is preserved under `analysis/WGCNA/archive/nonmanuscript_beta30/`. Its README states that it was not used in the submitted manuscript. No file in the current primary or reviewer-check directories depends on those module assignments.

## Remaining provenance limits

The exact upstream GeneID-to-symbol conversion script was not recovered. The exact processed FPKM matrix and checksum are therefore included. The beta-16 soft-threshold table was not archived in full; the release reports the manuscript-locked beta-16 R2 value without inventing unavailable slope or connectivity values.
