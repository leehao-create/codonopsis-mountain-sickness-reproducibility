# Differential expression and pathway analysis

The primary GSE103940 comparison is high altitude minus plain in 11 paired subjects using count-based TMM/voom-limma with subject in the design. The primary threshold is BH-adjusted P < 0.05 and absolute log2 fold change > 1.5. GSE75665 uses the subject-aware interaction `(AMS_high - AMS_plain) - (nonAMS_high - nonAMS_plain)` and pathway analysis of the full ranked statistic. GSE260910 is included only to reproduce the complete disease-batch confounding decision.

MSigDB v2025.1.Hs GMT files are license controlled and excluded. Place the Hallmark, Reactome and GO BP files named in the scripts in a directory and set `MSIGDB_DIR`.

The archived NCBI GeneID-to-symbol export is not redistributed because its exact upstream export procedure and snapshot date were not recovered. Run `scripts/download_ncbi_gene_mapping.sh` to create `data/reference/NCBI_GeneID_to_symbol.txt` from the current NCBI Homo sapiens `gene_info` file, or set `NCBI_GENE_MAP` to an authorized two-column tab-delimited mapping. The conversion uses NCBI `GeneID` and `Symbol` columns and preserves the header expected by the analysis scripts. A current download may differ from the archived mapping; its archived SHA-256 is retained in the top-level third-party input manifest for provenance, not as a claim of byte-identical reconstruction.

The GSE260910 audit optionally checks the public GEO family SOFT metadata. Download the family SOFT file from the [GSE260910 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE260910), save it as `GSE260910_family.soft.gz`, and set `GSE260910_SOFT` to that path. This script only reproduces the disease-batch confounding audit and does not perform HAPE differential expression.
