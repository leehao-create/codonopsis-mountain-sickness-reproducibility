# Differential expression and pathway analysis

The manuscript-locked GSE103940 result is reproduced from the archived GEO supplementary FPKM matrix by `scripts/01_GSE103940_manuscript_DEG.R`. The stored limma coefficient is `plain - high altitude`; signs are inverted when reporting the biological direction `high altitude - plain`. At BH-adjusted P < 0.05 and absolute log2 fold change > 1.5, this frozen table contains 2,782 genes: 118 increased and 2,664 decreased after high-altitude exposure. The archived workflow applied `log2(FPKM)`, replaced `log2(0)` values with zero, and used an unpaired condition-only limma design. These details are reported transparently because the repository is intended to reproduce the submitted numerical result, not to substitute the later 296-DEG count-based analysis.

The later paired count-based TMM/voom-limma and globin-sensitivity files are excluded from this module's public main-analysis release. Count/voom objects retained under `analysis/cell_composition/data/dependencies/` are scoped only to the separate reviewer-requested lineage-marker audit and must not be interpreted as the manuscript's primary GSE103940 DEG result.

GSE75665 uses the subject-aware interaction `(AMS_high - AMS_plain) - (nonAMS_high - nonAMS_plain)` and pathway analysis of the full ranked statistic. GSE260910 is the HAPE case-control dataset used in the manuscript analysis; the repository includes its post hoc audit documenting complete confounding between disease status and sequencing batch, a limitation that should be considered when interpreting HAPE-specific results.

MSigDB v2025.1.Hs GMT files are license controlled and excluded. Place the Hallmark, Reactome and GO BP files named in the scripts in a directory and set `MSIGDB_DIR`.

The archived NCBI GeneID-to-symbol export is not redistributed because its exact upstream export procedure and snapshot date were not recovered. Run `scripts/download_ncbi_gene_mapping.sh` to create `data/reference/NCBI_GeneID_to_symbol.txt` from the current NCBI Homo sapiens `gene_info` file, or set `NCBI_GENE_MAP` to an authorized two-column tab-delimited mapping. The conversion uses NCBI `GeneID` and `Symbol` columns and preserves the header expected by the analysis scripts. A current download may differ from the archived mapping; its archived SHA-256 is retained in the top-level third-party input manifest for provenance, not as a claim of byte-identical reconstruction.

The GSE260910 post hoc audit optionally checks the public GEO family SOFT metadata. Download the family SOFT file from the [GSE260910 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE260910), save it as `GSE260910_family.soft.gz`, and set `GSE260910_SOFT` to that path. This released script reproduces the disease-batch confounding audit; it does not replace or alter the HAPE case-control analysis reported in the manuscript.

## GSE103940 frozen files

- Input matrix: `data/processed/GSE103940_FPKM_gene_symbol_matrix.csv`
- Sample map: `data/metadata/GSE103940_subject_pairing.csv`
- Reproduction script: `scripts/01_GSE103940_manuscript_DEG.R`
- Frozen all-gene table: `results/GSE103940_manuscript/GSE103940_DEG_all_plain_minus_high.csv`
- Annotated workbook: `results/GSE103940_manuscript/GSE103940_DEG_annotated.xlsx`
- Exact current manuscript figure: `results/GSE103940_manuscript/Figure2_current_manuscript.pdf` and `.png`

The annotated workbook contains 2,665 rows including its header for high-altitude-down genes and 119 rows including its header for high-altitude-up genes, corresponding to 2,664 and 118 genes, respectively.

The exact R/limma package version used for the archived analysis was not recorded in the legacy files. The frozen all-gene table is therefore retained as the authoritative numerical output, and the reproduction script contains explicit count assertions to detect package-dependent drift.
