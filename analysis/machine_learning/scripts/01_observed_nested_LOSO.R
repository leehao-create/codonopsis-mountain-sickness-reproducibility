options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
module_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
suppressPackageStartupMessages({library(glmnet); library(randomForest); library(e1071)})
source(file.path(script_dir, "final_ml_utils.R"))

old <- readRDS(file.path(module_dir, "data", "dependencies",
                         "revised_ML_inputs_and_observed_fits.rds"))
out <- Sys.getenv("ML_OUTPUT_DIR", file.path(module_dir, "reproduced_results"))
dir.create(file.path(out, "GSE103940"), recursive = TRUE, showWarnings = FALSE)
x <- old$GSE103940$x; y <- old$GSE103940$y; subjects <- old$GSE103940$subjects
stopifnot(length(old$targets) == 54L, nrow(x) == 22L, ncol(x) == 39L,
          length(unique(subjects)) == 11L)

fit <- run_final_loso(x, y, subjects, 103940L, TRUE)
perf <- final_performance(fit$predictions)
paired <- paired_scores(fit$predictions)
paired_summary <- aggregate(correctly_ranked ~ method, paired,
                            function(z) c(n = sum(z), accuracy = mean(z)))
paired_summary <- data.frame(method = paired_summary$method,
                             subjects_correct = paired_summary$correctly_ranked[, "n"],
                             total_subjects = 11L,
                             paired_ranking_accuracy = paired_summary$correctly_ranked[, "accuracy"])

write.csv(fit$predictions, file.path(out, "GSE103940", "GSE103940_final_LOSO_OOF_predictions.csv"), row.names = FALSE)
write.csv(perf, file.path(out, "GSE103940", "GSE103940_final_OOF_performance.csv"), row.names = FALSE)
write.csv(fit$tuning, file.path(out, "GSE103940", "GSE103940_outer_fold_tuning.csv"), row.names = FALSE)
write.csv(fit$selections, file.path(out, "GSE103940", "GSE103940_outer_fold_feature_selections.csv"), row.names = FALSE)
write.csv(fit$leakage_audit, file.path(out, "GSE103940", "GSE103940_nested_validation_leakage_audit.csv"), row.names = FALSE)
write.csv(paired, file.path(out, "GSE103940", "GSE103940_paired_subject_score_differences.csv"), row.names = FALSE)
write.csv(paired_summary, file.path(out, "GSE103940", "GSE103940_paired_ranking_summary.csv"), row.names = FALSE)

boot <- subject_bootstrap(fit$predictions, 5000L, 4103940L)
ci <- do.call(rbind, lapply(unique(boot$method), function(method) {
  z <- boot[boot$method == method, ]
  data.frame(method = method, metric = "ROC_AUC",
             observed = perf$ROC_AUC[match(method, perf$method)],
             lower_95 = quantile(z$ROC_AUC, .025, na.rm = TRUE),
             median = quantile(z$ROC_AUC, .5, na.rm = TRUE),
             upper_95 = quantile(z$ROC_AUC, .975, na.rm = TRUE), row.names = NULL)
}))
write.csv(boot, file.path(out, "GSE103940", "GSE103940_subject_bootstrap_metrics.csv"), row.names = FALSE)
write.csv(ci, file.path(out, "GSE103940", "GSE103940_AUC_subject_bootstrap_95CI.csv"), row.names = FALSE)
saveRDS(list(input = list(targets = old$targets, x = x, y = y, subjects = subjects,
                          eligible_genes = old$GSE103940$eligible_genes),
             fit = fit, performance = perf, paired = paired,
             paired_summary = paired_summary, bootstrap_ci = ci),
        file.path(out, "final_ML_observed_fit.rds"))
capture.output(sessionInfo(), file = file.path(out, "final_ML_observed_sessionInfo.txt"))
print(perf); print(ci); print(paired_summary)
