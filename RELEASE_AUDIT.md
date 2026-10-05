# Public Release Audit

Audit date: 2026-10-05  
Scope: complete `reproducibility_repository` release tree. No biological or statistical analysis was rerun.

## Release decision

**Not ready for public release: one manuscript-method consistency blocker remains**

The repository files are internally aligned to the manuscript-locked GSE103940 numerical result (2,782 DEGs; 118 increased and 2,664 decreased at high altitude), and no file-integrity or redistribution blocker remains. However, `Revised2_Manuscript.docx` Methods 2.3 still describes a paired count-based TMM/voom-limma workflow, whereas the archived table yielding 2,782 DEGs used log2-transformed FPKM values and a condition-only limma model. Public release before correcting or explicitly explaining that mismatch would not satisfy the requested manuscript/repository consistency standard.

## Final audit matrix

| Area | Status | Finding |
|---|---|---|
| LICENSE | PASS | Standard MIT License is present and explicitly limited in the README to author-written repository code |
| Third-party licensing boundary | PASS | README and `DATA_PROVENANCE.md` state that third-party resources retain their original terms and are not relicensed |
| Human Cell Landscape objects | PASS | `blood1.rds` and `lung1.rds` were removed; accession, expected filenames, required object schema, archived sizes and SHA-256/MD5 values are retained |
| NCBI mapping | PASS | The untraceable archived mapping was removed; official source URL, conversion helper, expected output and archived digest are documented |
| STRING outputs | PASS | Three compact derived tables are retained with modification notice and attribution under STRING CC BY 4.0; query parameters and a current reconstruction helper are included |
| Other excluded inputs | PASS | Raw GeneCards, OMIM, TCMSP and HERB exports and MSigDB GMT files remain excluded and documented |
| External dependency documentation | PASS | GSE134355 objects, NCBI mapping, GSE260910 SOFT and MSigDB GMT acquisition/configuration are documented in top-level and module documentation |
| Manifest | PASS | `repository_file_manifest.csv` covers every tracked file except itself and the checksum list; every entry is final and marked releasable |
| SHA-256 checksums | PASS | `checksums.sha256` validates all 284 tracked files using repository-relative paths |
| Script syntax | PASS | 17 R scripts parsed; Python scripts compiled; shell scripts passed `bash -n` |
| Docking internal integrity | PASS | All 98 entries in the module-level docking checksum manifest validate after sanitizing non-scientific source-path remarks |
| CSV structure | PASS | All CSV files have consistent row widths |
| Absolute/local paths | PASS | No `/data3`, home-directory, workstation, volume or personal-account path remains |
| Credentials/privacy | PASS | No credential/private-key pattern or direct participant identifier was detected; author names in citation metadata and public accessions are intentional |
| Revised WGCNA only | PASS | No legacy-named WGCNA file is present; the repository contains only beta=30, 4,289 retained genes and three biological modules plus grey |
| GSE103940 release version | PASS (repository) | Main differential-expression files reproduce 2,782 DEGs (118 high-altitude increased; 2,664 decreased); the later 296-DEG analysis is excluded from the primary differential-expression module |
| GSE103940 manuscript-method agreement | BLOCKED | Manuscript Methods 2.3 describes paired count/TMM/voom, but the frozen 2,782-DEG table derives from log2(FPKM) and a condition-only limma design |
| Frozen ML | PASS | Final nested subject-level LOSO materials are retained; no new model or GSE260910 supervised analysis is included |
| Revised docking only | PASS | Final 2AA2/1GKC workflow and results are retained; no superseded docking conclusion was reintroduced |

## Third-party data resolution

### Removed from the public release

- `analysis/single_cell/data/processed/blood1.rds`
- `analysis/single_cell/data/processed/lung1.rds`
- `analysis/differential_expression/data/reference/NCBI_GeneID_to_symbol_archived.txt`

The source project copies were not altered. The public repository retains only provenance, archived checksums, expected object/file schemas, acquisition instructions and processing/downstream scripts. `.gitignore` rules prevent these user-supplied files from being committed accidentally.

### Retained with an explicit redistribution basis

- `analysis/network_pharmacology/data/derived/Fig1B_all_543_STRING_edges.csv`
- `analysis/network_pharmacology/data/derived/Fig1B_all_nodes_frozen_metrics.csv`
- `analysis/network_pharmacology/data/derived/Fig1C_top15_degree.csv`

STRING's official licensing page states that website, API and download outputs, including tables and data files, are available under Creative Commons Attribution 4.0. Repository documentation provides attribution, identifies the files as reformatted/derived, records the recoverable query settings, and notes that the exact historical STRING release was not recorded.

## Reproducibility boundaries

- The GSE103940 main-analysis script intentionally preserves the archived workflow that generated the manuscript-locked 2,782-DEG table. It does not claim that this workflow is paired or count based. The current manuscript Methods wording must be reconciled with this provenance before release.
- The exact R/limma version used for the legacy GSE103940 analysis was not recovered. The frozen table is retained as authoritative, and the script asserts the expected DEG counts so version drift fails visibly.
- Exact rerunning of the single-cell localization requires authorized custom `blood1.rds` and `lung1.rds` objects. GSE134355 source acquisition and the downstream object schema are documented, but the exact byte-identical source-to-RDS construction script was not recovered. This is disclosed rather than reconstructed by assumption.
- Current NCBI `gene_info` content may differ from the archived mapping. The provided helper creates a functional current mapping but does not claim byte identity.
- Current STRING results may differ from the frozen Figure 1 tables because the precise historical STRING release was not recorded. The versioned helper writes to a separate reconstruction directory and never overwrites frozen outputs.
- MSigDB files must be obtained by each user under the applicable MSigDB terms.
- A DOI and public repository URL can be added to `CITATION.cff` and `README.md` after upload. Their absence from a pre-upload tree is not a content-release blocker.

## Integrity commands

From the repository root:

```bash
sha256sum -c checksums.sha256
python tools/update_repository_manifest.py
```

The second command intentionally changes the manifest/checksum files and should only be used after an authorized repository edit; rerun the first command afterward.
