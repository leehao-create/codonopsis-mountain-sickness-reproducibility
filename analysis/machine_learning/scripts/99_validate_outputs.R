options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
module_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
out <- Sys.getenv("ML_VALIDATE_DIR", file.path(module_dir, "results"))
pred <- read.csv(file.path(out, "GSE103940", "GSE103940_final_LOSO_OOF_predictions.csv"))
audit <- read.csv(file.path(out, "GSE103940", "GSE103940_nested_validation_leakage_audit.csv"))
perf <- read.csv(file.path(out, "GSE103940", "GSE103940_final_OOF_performance.csv"))
perm <- read.csv(file.path(out, "GSE103940", "GSE103940_paired_permutation_null_metrics.csv"))
psum <- read.csv(file.path(out, "GSE103940", "GSE103940_paired_permutation_summary.csv"))
boot <- read.csv(file.path(out, "GSE103940", "GSE103940_subject_bootstrap_metrics.csv"))
stable <- read.csv(file.path(out, "ML_final_feature_stability.csv"))
paired <- read.csv(file.path(out, "GSE103940", "GSE103940_paired_subject_score_differences.csv"))

stopifnot(nrow(pred) == 22L, !anyDuplicated(pred$row_id),
          length(unique(pred$subject_id)) == 11L,
          all(table(pred$subject_id) == 2L), nrow(audit) == 11L,
          all(audit$train_test_subject_overlap == 0L),
          all(audit$test_pair_intact), !any(audit$inner_subject_split),
          setequal(perf$method, c("LASSO", "ElasticNet", "RF", "SVM", "Ensemble")),
          nrow(perm) == 499L * 5L,
          all(table(perm$method) == 499L),
          all(psum$empirical_permutation_P[psum$metric == "ROC_AUC"] < 0.05),
          nrow(boot) == 5000L * 5L,
          nrow(paired) == 11L * 5L, all(paired$correctly_ranked),
          nrow(stable) == 54L,
          identical(names(stable), c("gene", "LASSO_frequency",
                                     "ElasticNet_frequency", "RF_frequency",
                                     "SVM_frequency", "number_of_methods_supported")),
          all(stable$number_of_methods_supported[stable$gene %in% c("ATM", "SNCA", "BCL2")] == 4L),
          stable$number_of_methods_supported[stable$gene == "MMP9"] == 0L)

required_figures <- file.path(out, "fig5_components", paste0(
  c("Fig5A_final_ML_workflow", "Fig5B_GSE103940_OOF_ROC",
    "Fig5C_GSE103940_permutation_distributions",
    "Fig5D_GSE103940_feature_stability_heatmap",
    "Fig5E_GSE103940_consensus_candidates",
    "Fig5F_GSE103940_paired_score_differences"), ".pdf"))
stopifnot(all(file.exists(required_figures)),
          file.exists(file.path(out, "supplementary_GSE75665", "README.md")))

cat("PASS: 11 intact outer LOSO folds; zero subject overlap; no inner subject split\n")
cat("PASS: 22 unique OOF predictions and 11 intact paired comparisons\n")
cat("PASS: 499 paired permutations x 5 methods; all AUC P < 0.05\n")
cat("PASS: 5000 subject bootstrap iterations x 5 methods\n")
cat("PASS: 54-row final stability table with required columns\n")
cat("PASS: six Fig. 5 components and GSE75665 supplementary package\n")
