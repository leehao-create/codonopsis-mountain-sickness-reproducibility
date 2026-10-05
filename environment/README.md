# Software environment

The analyses were frozen under more than one R installation. Exact executed versions are retained in module-specific session files rather than represented by a synthetic lockfile:

- bulk differential expression and WGCNA: R 4.5.3;
- machine learning: R 4.5.3, glmnet 5.0, randomForest 4.7.1.2, e1071 1.7.17, ggplot2 4.0.2;
- single-cell localization: R 4.5.3, Matrix 1.7-5, ggplot2 4.0.2;
- cell-composition audit: R 4.6.1, edgeR 4.8.2, limma 3.66.0;
- docking: use `analysis/molecular_docking/environment.yml`; executed versions are in `software_versions.txt`.

Authoritative package snapshots are stored in the module-specific `sessionInfo.txt`, package-version CSV, and Conda environment files. No `renv.lock` was fabricated after the analyses. A clean-environment end-to-end rerun has not been performed as part of repository assembly.
