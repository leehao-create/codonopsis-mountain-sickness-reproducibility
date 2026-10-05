# Reproducibility Repository

This release-staging repository contains the frozen computational materials for the study **Integrative Multi-Omics and Network Pharmacology Analyses Reveal Potential Mechanisms of Codonopsis pilosula in Mountain Sickness**. It was assembled from the final analysis directories without rerunning any analysis.

## Study datasets

- GSE103940: high-altitude exposure transcriptomics (11 subjects, plain and high-altitude samples). The manuscript-locked DEG table reports 2,782 genes (118 increased and 2,664 decreased at high altitude); see the version note below.
- GSE75665: repeated-measures AMS study (10 subjects, plain and high-altitude samples; 5 AMS and 5 non-AMS).
- GSE260910: retained only for the documented disease-batch confounding audit; it is not used for inferential HAPE differential expression or supervised machine learning.
- GSE134355 / Human Cell Landscape reference objects: descriptive blood and lung cell-type localization only, not disease validation.

## Repository layout

```text
analysis/
  network_pharmacology/      frozen candidate, PPI and enrichment display data
  differential_expression/  paired/repeated-measures bulk analyses and GSEA
  WGCNA/                     final revised beta=30 signed network only
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
4. Run the revised WGCNA scripts in `analysis/WGCNA/scripts/`.
5. Run machine-learning scripts `01` through `04`; `99_validate_outputs.R` validates the frozen release outputs.
6. Run the single-cell localization only after supplying authorized `blood1.rds` and `lung1.rds` objects through `HCL_RDS_DIR`; run the cell-composition audit with the included GEO-derived inputs.
7. Create the docking environment from `analysis/molecular_docking/environment.yml`, activate it, and run `scripts/run_all.sh` only when full docking reproduction is desired.

No analysis was rerun during repository assembly. A syntax/static check does not constitute independent numerical reproduction.

## Authoritative GSE103940 manuscript version

The public main-analysis release follows the result frozen in the submitted manuscript: 2,782 DEGs, comprising 118 high-altitude-increased and 2,664 high-altitude-decreased genes. Its archived FPKM input, limma script, complete result table, annotated workbook and exact current Figure 2 are under `analysis/differential_expression/`. The later 296-DEG paired count/TMM/voom analysis is not released as the manuscript's primary differential-expression analysis. Count/voom dependencies that remain in the separate cell-composition and machine-learning modules serve those explicitly scoped analyses only.

## Authoritative WGCNA version

Only the final revised WGCNA is included: 4,289 retained genes, MAD threshold 0.50484, signed network and signed TOM, beta=30, `minModuleSize=20`, `deepSplit=2`, and `mergeCutHeight=0.25`. No candidate power reached R2 >= 0.85; beta=30 was selected by the documented fallback rule. The network contains three biological modules (turquoise, blue, brown) plus grey/unclassified genes. A global BH correction was applied across nine subject-aware module-effect tests, and no association remained significant. Legacy beta=16, 2,810-gene and 12-module outputs are excluded.

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

The repository tree has passed its local file-integrity and redistribution checks. The GSE103940 public main-analysis files now reproduce the manuscript-locked numerical result (2,782 DEGs), but the current manuscript Methods still describes the later paired count/TMM/voom workflow. This method/result wording mismatch is documented in `RELEASE_AUDIT.md` and should be resolved before public release. A persistent DOI/URL can be added to `CITATION.cff` and this README after upload.
