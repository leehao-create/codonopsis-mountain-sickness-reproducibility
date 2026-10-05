options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
module_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
suppressPackageStartupMessages({library(glmnet); library(randomForest); library(e1071)})
source(file.path(script_dir, "final_ml_utils.R"))

out <- Sys.getenv("ML_OUTPUT_DIR", file.path(module_dir, "reproduced_results"))
observed <- readRDS(file.path(out, "final_ML_observed_fit.rds"))
old_stability <- readRDS(file.path(module_dir, "data", "dependencies",
                                   "revised_ML_stability_results.rds"))
x <- observed$input$x; y <- observed$input$y; subjects <- observed$input$subjects
eligible <- observed$input$eligible_genes; targets <- observed$input$targets

rows <- vector("list", 600L)
k <- 0L
for (b in seq_len(200L)) {
  if (b == 1L || b %% 20L == 0L) message("Penalized stability ", b, "/200")
  set.seed(720260917L + b)
  sampled_subjects <- sample(sort(unique(subjects)), 9L, replace = FALSE)
  idx <- subjects %in% sampled_subjects
  lasso <- fit_penalized_nested(x[idx, , drop = FALSE], y[idx], subjects[idx],
                                1, "LASSO")
  elastic <- fit_penalized_nested(x[idx, , drop = FALSE], y[idx], subjects[idx],
                                  c(0.1, 0.3, 0.5, 0.7, 0.9), "Elastic Net")
  svm <- fit_svm_rfe_strict_nested(x[idx, , drop = FALSE], y[idx], subjects[idx],
                                   720260917L + b * 10000L)
  for (method in c("LASSO", "ElasticNet", "SVM")) {
    selected <- switch(method, LASSO = lasso$selected,
                       ElasticNet = elastic$selected, SVM = svm$selected)
    k <- k + 1L
    rows[[k]] <- data.frame(resample = b, method = method, gene = eligible,
                            selected = eligible %in% selected,
                            sampled_subjects = paste(sort(sampled_subjects), collapse = ";"))
  }
}
penalized_long <- do.call(rbind, rows)
write.csv(penalized_long, file.path(out, "GSE103940",
                                    "GSE103940_penalized_feature_stability_long.csv"),
          row.names = FALSE)

pen_freq <- aggregate(selected ~ method + gene, penalized_long, mean)
lasso_freq <- setNames(pen_freq$selected[pen_freq$method == "LASSO"],
                       pen_freq$gene[pen_freq$method == "LASSO"])
enet_freq <- setNames(pen_freq$selected[pen_freq$method == "ElasticNet"],
                      pen_freq$gene[pen_freq$method == "ElasticNet"])
svm_freq <- setNames(pen_freq$selected[pen_freq$method == "SVM"],
                     pen_freq$gene[pen_freq$method == "SVM"])
old <- old_stability$summary[old_stability$summary$dataset == "GSE103940", ]
rf_freq <- setNames(old$RF_frequency, old$gene)

final <- data.frame(
  gene = targets,
  LASSO_frequency = unname(lasso_freq[targets]),
  ElasticNet_frequency = unname(enet_freq[targets]),
  RF_frequency = unname(rf_freq[targets]),
  SVM_frequency = unname(svm_freq[targets]),
  stringsAsFactors = FALSE
)
final$number_of_methods_supported <- rowSums(
  cbind(final$LASSO_frequency >= 0.60, final$ElasticNet_frequency >= 0.60,
        final$RF_frequency >= 0.60, final$SVM_frequency >= 0.60), na.rm = TRUE)
stopifnot(nrow(final) == 54L)
write.csv(final, file.path(out, "ML_final_feature_stability.csv"), row.names = FALSE)
write.csv(final[order(-final$number_of_methods_supported,
                      -rowMeans(final[, 2:5], na.rm = TRUE)), ],
          file.path(out, "GSE103940", "GSE103940_feature_stability_ranked.csv"),
          row.names = FALSE)
saveRDS(list(penalized_long = penalized_long, final = final, threshold = 0.60,
             RF_source = "frozen Phase 5 200-resample results",
             SVM_source = "rerun with strict inner-fold feature ranking"),
        file.path(out, "final_ML_feature_stability.rds"))
capture.output(sessionInfo(), file = file.path(out, "final_ML_stability_sessionInfo.txt"))
print(final[final$number_of_methods_supported > 0L, ])
