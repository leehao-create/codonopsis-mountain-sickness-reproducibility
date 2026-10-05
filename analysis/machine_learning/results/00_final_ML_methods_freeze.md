# Final ML robustness methods freeze

Frozen before the final robustness models were run (2026-09-17).

## Scope

- Dataset: GSE103940 only for new model fitting.
- Candidate universe: the independently defined 54 Codonopsis pilosula-mountain
  sickness overlapping targets; no phenotype-dependent transcriptomic screening.
- Expression input: the frozen Phase 4 analysis-A TMM/voom matrix. Genes absent
  from that matrix are ineligible, leaving 39 predictors.
- Outcome: high altitude (1) versus plain/baseline (0), with 11 paired subjects.
- GSE75665 is not refitted. Its accepted Phase 5 results are packaged as
  supplementary exploratory evidence only.

## Validation and models

- Outer validation: leave one complete subject out (both paired samples).
- Inner grouping: samples from a training subject always share one inner fold.
- Preprocessing, zero-variance filtering, scaling, feature selection and tuning
  see outer-training subjects only.
- Penalized models use `glmnet` training-fold standardization. LASSO uses alpha=1.
  Elastic Net alpha is selected from {0.1, 0.3, 0.5, 0.7, 0.9}; lambda uses the
  one-standard-error rule after grouped inner CV.
- Random Forest retains the frozen Phase 5 specification. For linear SVM-RFE,
  each inner validation fold obtains scaling and its feature ranking from that
  inner fold's training subjects only; the selected cost/subset size is then
  refit on all outer-training subjects.
- The ensemble is the unweighted mean of LASSO, Elastic Net, RF and SVM
  probabilities. No ensemble weights are learned.

## Robustness analyses

- Paired permutation: 499 previously frozen, unique, non-observed within-subject
  label-swap assignments; empirical P=(1+#null >= observed)/(499+1).
- Subject bootstrap: 5,000 resamples of the 11 subjects with replacement, keeping
  each plain/high-altitude pair together. Percentile 95% CIs are reported.
- Feature stability: 200 resamples containing 9/11 subjects. Stable support is
  predefined as selection frequency >=0.60. LASSO and Elastic Net are fit in this
  phase; unchanged Phase 5 RF frequencies are carried forward, while SVM-RFE is
  rerun with the stricter inner-fold feature-ranking implementation.
- RF selection remains based on permutation importance exceeding both zero and
  its resample-specific 95th percentile null threshold; no Top-N rule is used.

No thresholds or tuning grids will be changed in response to the results.
