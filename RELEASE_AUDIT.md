# Public Release Audit

Audit date: 2026-10-13
Scope: complete `reproducibility_repository` release tree. No biological or statistical analysis was rerun.

## Release decision

**Ready for public release**

The repository files are internally aligned to the manuscript-locked GSE103940 numerical result (2,782 DEGs; 118 increased and 2,664 decreased at high altitude), and no file-integrity or redistribution blocker remains. The GSE103940 module README and provenance document transparently record the archived GEO supplementary FPKM, log2 transformation and condition-only limma workflow without claiming that all manuscript method wording and repository implementation are identical.

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
| SHA-256 checksums | PASS | `checksums.sha256` validates all 418 release files using repository-relative paths |
| Script syntax | PASS | 20 R scripts and 12 Python scripts parsed; shell scripts passed `bash -n` |
| Docking internal integrity | PASS | Frozen inputs, parameters, logs, redocking validation and candidate poses are present for NR3C2/2AA2, MMP9/1GKC, CA2/5E28 and BCL2/8HTS; CA2/BCL2 used Python-encoded parameters rather than standalone config text files |
| Compound/target provenance | PASS | Pair-level provenance for all reported docking pairs is retained in `analysis/molecular_docking/compound_target_provenance.csv`; the raw TCMSP mapping export remains excluded |
| Zenodo metadata | MANUAL CHECK | DOI `10.5281/zenodo.23228591` is recorded in README and CITATION.cff; Zenodo public accessibility requires manual browser confirmation |
| CSV structure | PASS | All CSV files have consistent row widths |
| Absolute/local paths | PASS | No absolute server, home-directory, workstation, volume or personal-account path remains |
| Credentials/privacy | PASS | No credential/private-key pattern or direct participant identifier was detected; author names in citation metadata and public accessions are intentional |
| Manuscript-primary WGCNA | PASS | Primary directory contains the 2,810-gene, beta-16, 12-module workflow; signed soft-threshold diagnostics are distinguished from the original module-construction call with signed TOM |
| WGCNA version separation | PASS | Reviewer checks are isolated; the 4,289-gene/beta-30 exploration is retained only under `archive/nonmanuscript_beta30/` and labeled as not used in the submitted manuscript |
| GSE103940 release version | PASS (repository) | Main differential-expression files reproduce 2,782 DEGs (118 high-altitude increased; 2,664 decreased); the later 296-DEG analysis is excluded from the primary differential-expression module |
| GSE103940 workflow provenance | PASS | The module README and version-provenance document record the archived GEO supplementary FPKM, gene filtering, log2 transformation, condition-only limma model, BH correction and frozen DEG threshold |
| Frozen ML | PASS | Final nested subject-level LOSO materials are retained; no new model or GSE260910 supervised analysis is included |
| Manuscript docking archive | PASS | Frozen 2AA2/1GKC and 5E28/8HTS workflows and results are retained; no superseded receptor structure or score was reintroduced |

## Reviewer 3 Minor Revision 26 material audit

| Reviewer requirement | Status | Exact repository location | Finding |
|---|---|---|---|
| Analysis scripts | PASS | `analysis/*/scripts/` | Author-written scripts are present for every released analysis module. |
| Processed intermediate files | PASS | Module-specific `data/processed/`, `data/dependencies/`, `prepared/` and `results/` directories | Redistribution-permitted matrices, fitted objects, prepared docking structures and frozen intermediate tables are included; restricted third-party files remain provenance-only. |
| Machine-learning random seeds | PASS | `analysis/machine_learning/README.md`; `results/00_final_ML_methods_freeze.md`; `results/supplementary/Table_S3C_model_specification.csv`; ML scripts | Deterministic folds and the recorded base/derived seeds for stability resampling, RF fitting, 499 permutations and 5,000 bootstraps are explicitly distinguished. |
| Frozen model/validation outputs | PASS | `analysis/machine_learning/results/GSE103940/`; `results/supplementary/`; `data/dependencies/` | OOF predictions, outer-fold tuning/selections, performance, permutation nulls, subject bootstrap metrics and frozen RDS objects are present. |
| Docking configurations | PASS WITH FORMAT NOTE | `analysis/molecular_docking/configs/`; `analysis/molecular_docking/ca2_bcl2/prepared/docking_boxes.csv`; CA2/BCL2 docking scripts and logs | NR3C2/MMP9 have standalone Vina configs. CA2/BCL2 passed the same fields through Python commands and did not create standalone config text files; no retrospective config was fabricated. |
| Docking poses/results/logs | PASS | `analysis/molecular_docking/results/` and `logs/`; `analysis/molecular_docking/ca2_bcl2/results/` and `logs/` | Candidate poses, redocking poses, all-pose RMSD tables, scores and execution logs are retained for 2AA2, 1GKC, 5E28 and 8HTS. |
| Compound/target provenance | PASS | `analysis/molecular_docking/compound_target_provenance.csv`; both docking `source_provenance.csv` files; `analysis/network_pharmacology/` | Pair-level mapping source, receptor identifiers, ligand identifiers/CIDs and retained input locations are documented without redistributing the raw TCMSP export. |
| Repository manifest | PASS | `repository_file_manifest.csv` | All public-release files except the manifest/checksum files themselves are indexed. |
| SHA-256 checksums | PASS | `checksums.sha256`; module-level input checksum lists | Repository-wide and docking-input integrity records validate. |
| Zenodo DOI metadata | MANUAL CHECK | `README.md`; `CITATION.cff`; this audit | DOI `10.5281/zenodo.23228591` is recorded. Zenodo public accessibility requires manual browser confirmation. |

The CA2/BCL2 archive contains no standalone Vina `.txt` config because none was produced by the original Python-driven workflow. It also contains no separate MMFF94 minimization log; the recovered CA2/BCL2 scripts/logs document PubChem 3D input, pH 7.4 protonation and Gasteiger/PDBQT preparation. These absent records were not recreated or inferred.

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

- The GSE103940 main-analysis script intentionally preserves the archived workflow that generated the manuscript-locked 2,782-DEG table. It does not claim that this workflow is paired or count based, nor that all manuscript method wording and repository implementation are identical.
- The WGCNA primary script intentionally preserves the original module-construction call. Signed soft-threshold diagnostics and a signed TOM do not establish that `networkType` was explicitly set for `blockwiseModules`; repository documentation states this boundary directly.
- The exact R/limma version used for the legacy GSE103940 analysis was not recovered. The frozen table is retained as authoritative, and the script asserts the expected DEG counts so version drift fails visibly.
- Exact rerunning of the single-cell localization requires authorized custom `blood1.rds` and `lung1.rds` objects. GSE134355 source acquisition and the downstream object schema are documented, but the exact byte-identical source-to-RDS construction script was not recovered. This is disclosed rather than reconstructed by assumption.
- Current NCBI `gene_info` content may differ from the archived mapping. The provided helper creates a functional current mapping but does not claim byte identity.
- Current STRING results may differ from the frozen Figure 1 tables because the precise historical STRING release was not recorded. The versioned helper writes to a separate reconstruction directory and never overwrites frozen outputs.
- MSigDB files must be obtained by each user under the applicable MSigDB terms.
- The archived release DOI is `10.5281/zenodo.23228591` and the GitHub repository URL is recorded in both release documentation and citation metadata. Zenodo public accessibility requires manual browser confirmation.

## Integrity commands

From the repository root:

```bash
sha256sum -c checksums.sha256
python tools/update_repository_manifest.py
```

The second command intentionally changes the manifest/checksum files and should only be used after an authorized repository edit; rerun the first command afterward.
