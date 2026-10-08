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
| Gene handling | supplementary files merged by gene symbol; genes with zero expression across all samples removed |
| Frozen transformation | log2(FPKM), with log2(0) replaced by zero |
| Frozen model | condition-only limma model |
| Multiple-testing correction | Benjamini-Hochberg |
| Exact manuscript figure | `analysis/differential_expression/results/GSE103940_manuscript/Figure2_current_manuscript.*` |

The later 296-DEG count-based paired TMM/voom-limma analysis is not the manuscript's primary GSE103940 analysis and is excluded from the public differential-expression module. Count/voom objects retained for cell-composition and machine-learning analyses remain confined to those modules and are labelled accordingly.

GSE103940 differential-expression analysis used the GEO supplementary FPKM expression files. Expression data were merged by gene symbol, genes with zero expression across all samples were removed, and FPKM values were log2-transformed, with non-finite values replaced by zero. Differential expression between baseline and high-altitude exposure samples was assessed using limma with condition as the model term. Benjamini-Hochberg correction was applied, and adjusted P < 0.05 together with |log2FC| > 1.5 was used to define differentially expressed genes.

This provenance statement documents the archived workflow without changing the manuscript, code, frozen results or Figure 2. It does not claim that all manuscript method wording and repository implementation are identical.
