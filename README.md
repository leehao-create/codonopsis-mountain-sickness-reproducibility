# Reproducibility Repository

This release-staging repository contains the frozen computational materials for the study **Integrative Multi-Omics and Network Pharmacology Analyses Reveal Potential Mechanisms of Codonopsis pilosula in Mountain Sickness**. It was assembled from the final analysis directories without rerunning any analysis.

## Study datasets

- GSE103940: high-altitude exposure transcriptomics (11 subjects, plain and high-altitude samples). The manuscript-locked DEG table reports 2,782 genes (118 increased and 2,664 decreased at high altitude); see the version note below.
- GSE75665: repeated-measures AMS study (10 subjects, plain and high-altitude samples; 5 AMS and 5 non-AMS).
- GSE260910: HAPE case-control transcriptomic dataset used in the manuscript analysis. A post hoc audit identified complete confounding between disease status and sequencing batch; this limitation is documented in the repository and should be considered when interpreting HAPE-specific results.
- GSE134355 / Human Cell Landscape reference objects: descriptive blood and lung cell-type localization only, not disease validation.

## Repository layout

```text
analysis/
  network_pharmacology/      frozen candidate, PPI and enrichment display data
  differential_expression/  paired/repeated-measures bulk analyses and GSEA
  WGCNA/                     manuscript-primary workflow plus isolated reviewer checks and archive
  machine_learning/          frozen nested subject-level LOSO analysis
  single_cell/               reference-atlas localization
  cell_composition/          lineage-marker expression audit
  molecular_docking/         final validated 2AA2/1GKC redocking workflow
environment/                 recorded software versions and setup notes
external_data/               manifests for excluded third-party inputs
tools/                       manifest/checksum maintenance
```

Each module contains `scripts/`, inputs or frozen dependencies where redistribution is currently considered permissible, and final `results/`. Scripts default to a module-local `reproduced_results/` directory so they do not overwrite the frozen results.

## Reproduction order

1. Verify the release with `sha256sum -c checksums.sha256`.
2. Read `DATA_PROVENANCE.md` and `external_data/THIRD_PARTY_INPUTS.csv`; obtain excluded third-party inputs, especially MSigDB v2025.1 gene-set files and the authorized Human Cell Landscape objects.
3. Run differential-expression scripts in numerical order.
4. Run or inspect the manuscript-primary WGCNA materials in `analysis/WGCNA/manuscript_primary/`; reviewer-only diagnostics are isolated under `analysis/WGCNA/reviewer_checks/`.
5. Run machine-learning scripts `01` through `04`; `99_validate_outputs.R` validates the frozen release outputs.
6. Run the single-cell localization only after supplying authorized `blood1.rds` and `lung1.rds` objects through `HCL_RDS_DIR`; run the cell-composition audit with the included GEO-derived inputs.
7. Create the docking environment from `analysis/molecular_docking/environment.yml`, activate it, and run `scripts/run_all.sh` only when full docking reproduction is desired.

No analysis was rerun during repository assembly. A syntax/static check does not constitute independent numerical reproduction.

## Authoritative GSE103940 manuscript version

The public main-analysis release follows the result frozen in the submitted manuscript: 2,782 DEGs, comprising 118 high-altitude-increased and 2,664 high-altitude-decreased genes. Its archived FPKM input, limma script, complete result table, annotated workbook and exact current Figure 2 are under `analysis/differential_expression/`. The later 296-DEG paired count/TMM/voom analysis is not released as the manuscript's primary differential-expression analysis. Count/voom dependencies that remain in the separate cell-composition and machine-learning modules serve those explicitly scoped analyses only.

## Manuscript-primary WGCNA version

The manuscript-primary GSE75665 WGCNA used 2,810 genes after `MAD > 1.5`, beta = 16, and produced 12 modules including grey. Soft-threshold diagnostics were evaluated using a signed-network setting, while module construction used the original `blockwiseModules` workflow with a signed TOM. The original module-construction call did not explicitly set `networkType` and is preserved without altering its behavior. The signed soft-threshold fit at beta 16 was R2 = 0.0419 and did not satisfy the conventional R2 >= 0.85 criterion.

The original 12 x 4 group-indicator module-trait matrices and reviewer-requested binary-AMS check are reported separately. No module association remained significant after the applicable BH correction; module-size-aware target enrichment was also non-significant. MMP9 is retained only as exploratory module context and is not presented as a high-connectivity hub. The separate 4,289-gene/beta-30 reanalysis is retained under `analysis/WGCNA/archive/nonmanuscript_beta30/` and is explicitly not used in the submitted manuscript.

## Analysis modules

1. Network pharmacology
2. Differential expression and enrichment
3. WGCNA manuscript-primary analysis
4. Machine-learning nested LOSO validation
5. Single-cell reference analysis
6. Cell-composition reviewer analysis
7. Molecular docking
8. Reviewer-specific sensitivity checks

Results under `reviewer_checks/` are diagnostic or sensitivity analyses and are not manuscript-primary findings.

## Integrity and release metadata

- `repository_file_manifest.csv` describes every tracked file, its role, size and SHA-256 digest.
- `checksums.sha256` contains repository-relative checksums and intentionally excludes itself and the CSV manifest.
- `CITATION.cff` supplies repository citation metadata.
- `RELEASE_AUDIT.md` records the final pre-release decision and unresolved author/licensing actions.

Regenerate integrity files after an authorized change:

```bash
python tools/update_repository_manifest.py
sha256sum -c checksums.sha256
```

## License and third-party data

The MIT License in `LICENSE` applies only to analysis, validation and rendering code authored for this repository. It does not relicense data, database content, structures, annotations or other material obtained from third parties.

GEO, Human Cell Landscape, STRING, NCBI, MSigDB, RCSB PDB, PubChem, GeneCards, OMIM, TCMSP, HERB and all other third-party resources remain governed by their original licenses, access conditions and terms of use. Users are responsible for obtaining those materials from the cited source and complying with the applicable terms. This repository makes no grant of rights over third-party content.

Raw GeneCards, OMIM, TCMSP and HERB exports, MSigDB GMT files, custom Human Cell Landscape RDS objects and the archived NCBI GeneID mapping are excluded. Source/accession, retrieval instructions, expected filenames and available archived checksums are recorded in `DATA_PROVENANCE.md`, `external_data/THIRD_PARTY_INPUTS.csv`, and module READMEs. Three compact STRING-derived tables are retained with attribution because STRING explicitly distributes website/API/download outputs under CC BY 4.0; their historical release was not recorded and the limitation is documented.

## Release status

A persistent DOI/URL can be added to `CITATION.cff` and this README after upload.
