# GSE103940 Version Alignment

This release follows the GSE103940 result currently reported in `Revised2_Manuscript.docx`.

| Item | Authoritative release value |
|---|---|
| Reported biological contrast | high altitude minus plain |
| Stored result-table contrast | plain minus high altitude |
| Significance rule | BH-adjusted P < 0.05 and absolute log2 fold change > 1.5 |
| Total DEGs | 2,782 |
| Increased at high altitude | 118 |
| Decreased at high altitude | 2,664 |
| Input data | archived GEO supplementary FPKM values |
| Frozen transformation | log2(FPKM), with log2(0) replaced by zero |
| Frozen model | condition-only limma model |
| Exact manuscript figure | `analysis/differential_expression/results/GSE103940_manuscript/Figure2_current_manuscript.*` |

The later 296-DEG count-based paired TMM/voom-limma analysis is not the manuscript's primary GSE103940 analysis and is excluded from the public differential-expression module. Count/voom objects retained for cell-composition and machine-learning analyses remain confined to those modules and are labelled accordingly.

The version alignment is numerical and provenance-based. It does not change the manuscript or retrospectively change the archived statistical workflow.
