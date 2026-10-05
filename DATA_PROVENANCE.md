# Data Provenance and Redistribution

| Content | Source | Repository treatment | Release status |
|---|---|---|---|
| GSE103940 processed FPKM and sample map | NCBI GEO GSE103940 | Gene-symbol FPKM matrix used for the manuscript-locked 2,782-DEG result and author-curated pairing map included | Included; cite the GEO accession |
| GSE75665 counts, FPKM and sample map | NCBI GEO GSE75665 | Included as accession-derived inputs and author-curated repeated-measures map | Included; cite the GEO accession |
| GSE260910 metadata | NCBI GEO GSE260910 | Sample map and confounding result included; no inferential HAPE result included | Releasable as derived audit material |
| MSigDB Hallmark, Reactome and GO BP GMT | MSigDB v2025.1.Hs | Not included | Users must obtain the named GMT files under the applicable MSigDB terms and set `MSIGDB_DIR` |
| 54 overlapping targets | Derived from archived network-pharmacology workflow | Compact derived gene list included; raw database exports excluded | Author should cite all contributing databases |
| GeneCards and OMIM exports | GeneCards / OMIM | Excluded | Do not redistribute raw exports |
| TCMSP and HERB exports | TCMSP / HERB | Excluded; derived compound table retained | Do not redistribute raw exports without permission |
| STRING/Cytoscape display tables | Frozen STRING-derived PPI results | Three compact Figure 1 plot tables retained with attribution | Included under STRING CC BY 4.0; exact historical STRING release was not recorded |
| Blood/lung RDS objects | Human Cell Landscape reference atlas, associated with GSE134355 | Custom processed objects removed; frozen derived localization tables and object checksums retained | Not redistributed; obtain source matrices from GSE134355 and provide authorized objects through `HCL_RDS_DIR` |
| NCBI GeneID-to-symbol mapping | NCBI Gene `Homo_sapiens.gene_info.gz` | Archived mapping removed; current download/conversion helper and archived digest retained | Not redistributed; run `analysis/differential_expression/scripts/download_ncbi_gene_mapping.sh` or set `NCBI_GENE_MAP` |
| GSE260910 family SOFT | NCBI GEO GSE260910 | Not included; sample map and final confounding audit are retained | Download from GEO and set `GSE260910_SOFT`; used only for metadata verification |
| PDB structures 2AA2 and 1GKC | RCSB PDB | Included with identifiers and provenance | Preserve RCSB/PDB attribution |
| Candidate ligand SDF files | Archived PubChem 3D records; CIDs in `source_provenance.csv` | Included | Preserve PubChem attribution |

## External acquisition and reconstruction

### Human Cell Landscape

- Source: NCBI GEO GSE134355, <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134355>.
- Source content: per-sample raw UMI digital-expression matrices for Human Cell Landscape tissues; use Adult Peripheral Blood and Adult/Fetal Lung records for the two released localization datasets.
- Expected downstream filenames: `blood1.rds` and `lung1.rds`, supplied in a user-controlled directory through `HCL_RDS_DIR`.
- Required object schema and archived checksums: `analysis/single_cell/README.md` and `analysis/single_cell/external_data_manifest.csv`.
- Limitation: the exact byte-identical source-to-RDS construction script was not recovered. No replacement thresholds or annotations were inferred. Frozen localization output tables are included, but the excluded RDS objects are not sublicensed by this repository.

### STRING

- Source and license: STRING v12.5 access/licensing page, <https://version-12-5.string-db.org/cgi/access?footer_active_subpage=licensing>; STRING states that tables and data files from its website, APIs and downloads are available under CC BY 4.0.
- Recoverable query: 54 supplied human genes; species `9606`; functional network; `required_score=400`; no added interactors.
- Historical limitation: the precise STRING release used for the frozen tables was not archived. The included outputs contain 54 nodes, 543 edges and a minimum combined score of 0.400.
- Reconstruction helper: `analysis/network_pharmacology/scripts/retrieve_current_STRING_network.py` uses the versioned v12.5 API and writes current outputs under `reconstructed_data/`; current results may differ from the frozen publication tables.
- Retained files: `Fig1B_all_543_STRING_edges.csv`, `Fig1B_all_nodes_frozen_metrics.csv`, and `Fig1C_top15_degree.csv`. These are modified/reformatted STRING-derived outputs and require STRING attribution.

### NCBI Gene mapping

- Source: <https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz>.
- Query/scope: Homo sapiens Gene records; the helper extracts `GeneID` and `Symbol`.
- Expected output: `analysis/differential_expression/data/reference/NCBI_GeneID_to_symbol.txt`, or any equivalent two-column file supplied via `NCBI_GENE_MAP`.
- Helper: `analysis/differential_expression/scripts/download_ncbi_gene_mapping.sh`.
- Archived mapping digest: SHA-256 `510407a793c744f7929122ca55e8ff57d6530d2163d44d820d9e445acba19207`. The snapshot date and exact original export procedure were not recovered, so a current mapping is not claimed to be byte-identical.

### GSE260910 metadata

- Source: NCBI GEO GSE260910, <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE260910>.
- Input: GEO family SOFT, expected filename `GSE260910_family.soft.gz`.
- Configuration: set `GSE260910_SOFT` to the downloaded file before running `analysis/differential_expression/scripts/02_GSE260910_batch_audit.R`.
- Scope: metadata and design-rank audit only; the repository does not perform inferential HAPE differential expression or supervised machine learning with this dataset.

No controlled-access clinical data or direct participant identifiers are included. Public GSM accessions and study pseudonyms are retained where required to reproduce paired or repeated-measures designs.

### GSE103940 version boundary

The main differential-expression release uses the archived GEO supplementary FPKM values in `analysis/differential_expression/data/processed/GSE103940_FPKM_gene_symbol_matrix.csv`. The reproduction script preserves the archived transformation and condition-only limma design and asserts the manuscript counts of 2,782 total DEGs, 118 increased and 2,664 decreased at high altitude. A later paired count/TMM/voom analysis yielding 296 DEGs is excluded from the main differential-expression release. Frozen count/voom dependencies retained for the reviewer-requested cell-composition audit are explicitly scoped to that audit and are not a substitute primary DEG result.

## Licensing boundary

The repository MIT License applies only to author-written code. All third-party data and resources remain subject to their source terms. Inclusion of a derived output does not transfer ownership or relicense its upstream source; the repository provides attribution and modification notices where applicable.
