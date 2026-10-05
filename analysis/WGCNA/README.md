# Final Revised WGCNA: GSE75665

## Data source and design

- GEO accession: **GSE75665**.
- Expression input: archived FPKM matrix, 24,901 genes x 20 samples.
- Subjects: 10, with paired plain and high-altitude samples for each subject.
- AMS status: 5 AMS and 5 non-AMS subjects.
- Included processed matrix: `data/processed/GSE75665_WGCNA_input_FPKM.csv`.
- Sample map: `data/metadata/GSE75665_subject_pairing.csv`.
- Predefined candidate-target list: `data/reference/Codonopsis_54_overlap_targets.csv`.

The included FPKM matrix is the exact matrix used by the final revised WGCNA and is derived from the public GEO processed expression resource for GSE75665. The archived upstream GEO file was `GSE75665_norm_counts_FPKM_GRCh38.p13_NCBI.tsv.gz`; the final input uses gene symbols as row identifiers and the 20 GSM accessions as columns. A separate, fully automated GeneID-to-symbol conversion script was not recovered from the archive; therefore, the exact final matrix and its checksum are provided rather than reconstructing undocumented conversion steps.

## Preprocessing and gene filtering

The authoritative analysis script is `scripts/03_GSE75665_signed_WGCNA.R`.

1. FPKM values were transformed as `log2(FPKM + 1)`.
2. Genes were expression-eligible when FPKM >= 1 in at least 5 of 20 samples.
3. Zero-MAD genes were removed.
4. The primary phenotype-blind variance filter retained the upper 50% of nonzero MAD values.
5. The observed median MAD threshold was **0.50484**.
6. **4,289 genes** were retained.
7. Sample outliers required both connectivity Z < -2.5 and robust distance Z > 3.5. No sample met this joint criterion; all 20 samples were retained.

## Network parameters

| Parameter | Final value |
|---|---|
| Network type | signed |
| TOM type | signed |
| Correlation | Pearson |
| Soft-threshold power | beta = 30 |
| Scale-free R2 at beta 30 | 0.103228 |
| Slope at beta 30 | -0.178930 |
| Mean connectivity at beta 30 | 753.020642 |
| Prespecified R2 threshold | 0.85, not reached |
| Fallback rule | Maximum R2 among candidate powers with negative slope and mean connectivity >= 1 |
| minModuleSize | 20 |
| deepSplit | 2 |
| mergeCutHeight | 0.25 |
| Random seed | 75665 |

Powers 1-30 were the candidate range. Powers 31-50 were evaluated only as a diagnostic extension; no power reached R2 >= 0.85. Beta 30 was selected by the prespecified deterministic fallback rule and **did not** meet the scale-free topology threshold.

## Module-trait analysis

The final network contained turquoise (4,085 genes), blue (94 genes), and brown (52 genes) modules, plus 58 grey/unclassified genes.

Each module eigengene was analyzed using the subject-aware mixed model:

```text
ME ~ AMS_status + altitude_status + AMS_status:altitude_status + (1 | subject_id)
```

The three module eigengenes and three fixed effects produced nine module-effect tests. A single global Benjamini-Hochberg correction was applied across all nine nominal P values. No association remained significant at BH-adjusted P < 0.05; the minimum global adjusted P was 0.354750.

Module-size-aware enrichment of the 54 predefined *Codonopsis pilosula* targets was assessed with one-sided Fisher exact tests and BH correction across the three biological modules. No module was significantly enriched after correction. Target assignments are exploratory co-expression context only.

MMP9 was assigned to the turquoise module (`kME = 0.723949`, `GS = 1.743908`, `kWithin = 34.078389`) and did not satisfy the prespecified hub criteria. Its module assignment is exploratory only.

## Reproduction order

1. Confirm that the three input files under `data/` match `repository_file_manifest.csv` and `checksums.sha256`.
2. Install the R/package versions recorded in `results/WGCNA_sessionInfo.txt` (R 4.5.3, WGCNA 1.74, limma 3.66.0, lmerTest 3.2-1, lme4 2.0-6, Matrix 1.7-5, and ggplot2 4.0.2).
3. Run `scripts/03_GSE75665_signed_WGCNA.R` to recreate the WGCNA results in `reproduced_results/`. Set `WGCNA_OUTPUT_DIR` only if another output location is required.
4. Run `scripts/01_build_revised_WGCNA_figure_and_tables.R` to render Figure 3 and CSV mirrors from the frozen `results/` directory. Set `WGCNA_RESULT_DIR` or `WGCNA_FIGURE_DIR` only when required.

The public repository copies differ from the archived originals only in filesystem plumbing: local absolute paths were replaced with repository-relative paths, outputs were redirected away from frozen results, and the unrelated merge into a master Supplementary workbook was removed. Analytical parameters and calculations were not changed.

## Included results

The `results/` directory contains one authoritative copy of the final configuration, filter and QC summaries, power diagnostics and decision, module assignments and sizes, eigengenes, subject-aware module-trait results, gene-level MM/GS/connectivity metrics, candidate metrics, formal target enrichment, hub table, saved final network object, outlier decision, and software versions.

S4A-S4H CSV mirrors were not copied because they duplicate these authoritative result files. Sensitivity networks involving beta = 16 and all legacy 2,810-gene/12-module outputs were excluded.
