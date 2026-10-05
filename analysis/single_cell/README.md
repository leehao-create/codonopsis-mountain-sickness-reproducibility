# Single-cell reference localization

The analysis reports average normalized expression and percent expressed by annotated cell type for predefined candidate genes. It uses blood and lung reference objects associated with the Human Cell Landscape accession GSE134355, does not recluster cells, and must not be interpreted as disease validation.

## External objects

`blood1.rds` and `lung1.rds` are custom processed objects and are **not redistributed**. Their archived byte counts and SHA-256/MD5 digests are recorded in `external_data_manifest.csv` and the frozen `results/single_cell_input_md5.csv`. Obtain the source digital-expression matrices from the [GSE134355 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134355). GEO identifies the study as Human Cell Landscape single-cell expression profiling and provides raw UMI digital-expression matrices per sample.

The exact byte-identical source-to-RDS construction script was not recovered, so the archived RDS objects cannot be reconstructed exactly from this release alone. The recovered workflow used Seurat normalization, variable-feature selection, scaling, PCA and Harmony integration, followed by manual cell-type annotation; the final downstream script requires each object to contain an RNA assay with `counts` and normalized `data`, and metadata fields `celltype`, `donor_id`, `batch`, `disease`, `assay`, and optionally `tissue_original`. This limitation is reported rather than filled with inferred processing steps.

To rerun the released localization step, place the two objects under any directory with the expected names and set:

```bash
export HCL_RDS_DIR=/path/to/authorized/hcl_objects
Rscript scripts/01_single_cell_candidate_localization.R
```

The expected outputs are the candidate-expression, coverage, major-cell-type, dataset-summary and donor-count tables plus the candidate dot plot in `reproduced_results/`. Frozen derived results remain included for inspection without redistributing the RDS objects.
