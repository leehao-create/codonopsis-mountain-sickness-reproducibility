# GSE75665 WGCNA materials

## Version boundary

`manuscript_primary/` is the only manuscript-primary WGCNA version. It contains the original 2,810-gene, beta-16 workflow that produced the 12-module results used in the submitted manuscript.

Soft-threshold diagnostics were evaluated using a signed-network setting, while module construction used the original `blockwiseModules` workflow with a signed TOM. The original module-construction call did not explicitly set `networkType`; this repository preserves that call and does not describe it as an explicitly signed adjacency analysis.

`reviewer_checks/` contains later binary-AMS, multiplicity, target-enrichment and MMP9 diagnostic checks. These analyses do not redefine the primary modules.

`archive/nonmanuscript_beta30/` preserves the 4,289-gene, beta-30 exploratory reanalysis for provenance. It was not used in the submitted manuscript and is neither current nor authoritative.

## Data source

- GEO accession: GSE75665.
- Processed input: `data/processed/GSE75665_WGCNA_input_FPKM.csv` (24,901 genes x 20 samples).
- Subjects: 10, each with paired plain and high-altitude samples; 5 AMS and 5 non-AMS subjects.
- Metadata: `data/metadata/GSE75665_subject_pairing.csv`.
- Predefined target reference: `data/reference/Codonopsis_54_overlap_targets.csv`.

The exact processed FPKM matrix is included because the archived GeneID-to-symbol conversion was not fully scripted. The matrix is derived from the public GEO processed resource `GSE75665_norm_counts_FPKM_GRCh38.p13_NCBI.tsv.gz`.

## Manuscript-primary parameters

| Parameter | Value |
|---|---|
| Samples | 20 |
| Genes after missing/zero-total filtering | 20,751 |
| MAD filter | MAD > 1.5 |
| Retained genes | 2,810 |
| Soft-threshold diagnostic setting | signed |
| Soft-threshold power | beta = 16 |
| Scale-free topology fit at beta 16 | R2 = 0.0419 |
| Conventional criterion | R2 >= 0.85, not reached |
| Module construction | original `blockwiseModules` call; `networkType` not explicitly set |
| TOM type | signed |
| Correlation | Pearson |
| Minimum module size | 20 |
| Merge cut height | 0.05 |
| Reassign threshold | 0 |
| Maximum block size | 2,810 |
| Modules, including grey | 12 |

## Statistical interpretation

The original module-trait output is a 12-module by four-group-indicator analysis (`high_AMS`, `high_non-AMS`, `low_AMS`, and `low_non-AMS`). Global BH correction across those 48 values gives a minimum adjusted P value of 0.387. A reviewer-requested binary AMS check tested 12 module eigengenes and gives a minimum BH-adjusted P value of 0.0607. Neither analysis supports an FDR-significant module association.

Module-size-aware candidate-target enrichment was non-significant after BH correction. MMP9 was assigned to turquoise and had `kME = 0.852`, traditional binary-trait `GS = 0.4423` (`P = 0.0508`), signed-adjacency diagnostic `kWithin = 226.2`, and within-module rank 967/1,147. These post hoc diagnostics are exploratory; MMP9 is not interpreted as a high-connectivity hub.

## Reproduction and verification

The byte-preserved archived code is `manuscript_primary/scripts/00_original_WGCNA_code_as_archived.R.txt`. It contains two pasted console-output lines and is retained for provenance rather than execution. The repository-portable transcription is `manuscript_primary/scripts/01_GSE75665_manuscript_primary_WGCNA.R`; it starts from the included processed matrix and intentionally preserves the original analytical calls, including the implicit `blockwiseModules` network-type default.

`manuscript_primary/scripts/02_format_frozen_WGCNA_outputs.R` formats archived outputs into public CSV tables and recomputes only BH columns from already saved P values; it does not rebuild a network. From the repository root:

```bash
Rscript analysis/WGCNA/manuscript_primary/scripts/02_format_frozen_WGCNA_outputs.R .
```

The reviewer-check plotting script writes only to `reviewer_checks/reproduced_figures/` and does not overwrite frozen result tables.
