# Machine-learning candidate prioritization

This is the frozen GSE103940 candidate-prioritization analysis, not a diagnostic model. It uses 39 expression-eligible genes from 54 predefined targets, outer leave-one-subject-out validation, grouped inner resampling, LASSO, Elastic Net, random forest and linear SVM-RFE. Scaling, tuning and feature selection occur inside outer-training subjects. Frozen OOF predictions, bootstrap confidence intervals, paired-label permutation results, stability results and Figure 5 plot data are included.

Run scripts `01` through `04` in order to write a separate `reproduced_results/` tree. Script `99_validate_outputs.R` checks the frozen `results/` tables. The carried-forward RF stability source and permutation codes are retained in `data/dependencies/` and documented in `00_final_ML_methods_freeze.md`.
