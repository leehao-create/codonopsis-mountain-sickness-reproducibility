# Machine-learning candidate prioritization

This is the frozen GSE103940 candidate-prioritization analysis, not a diagnostic model. It uses 39 expression-eligible genes from 54 predefined targets, outer leave-one-subject-out validation, grouped inner resampling, LASSO, Elastic Net, random forest and linear SVM-RFE. Scaling, tuning and feature selection occur inside outer-training subjects. Frozen OOF predictions, bootstrap confidence intervals, paired-label permutation results, stability results and Figure 5 plot data are included.

Run scripts `01` through `04` in order to write a separate `reproduced_results/` tree. Script `99_validate_outputs.R` checks the frozen `results/` tables. The carried-forward RF stability source and permutation codes are retained in `data/dependencies/` and documented in `00_final_ML_methods_freeze.md`.

## Validation and random-seed rules

- Outer validation is deterministic 11-fold subject-level LOSO; each held-out fold contains both samples from one subject.
- Grouped fivefold inner partitions are deterministic functions of the sorted training-subject identifiers.
- Scaling, feature selection and hyperparameter tuning are performed using outer-training data only.
- LASSO, Elastic Net and SVM-RFE stability use 200 nine-subject resamples with replicate-specific seeds derived from `720260917`.
- Random-forest stability uses base seed `4103940`; its resample and observed/null model seeds are derived in `upstream_ml_utils.R`.
- The 499 unique nonidentity paired-label assignments were sampled with seed `5103940`; model fitting under each assignment uses deterministic seeds derived from base seed `103940`.
- The 5,000 subject-level bootstrap resamples use seed `4103940` and resample frozen outer-fold predictions without retraining models.
- Random-forest importance is out-of-bag permutation importance from `randomForest::importance(type = 1, scale = FALSE)`, i.e. unscaled mean decrease in accuracy.

Exact tuning grids and derived-seed formulas are recorded in `results/supplementary/Table_S3C_model_specification.csv` and the scripts. No model was retrained during repository assembly.
