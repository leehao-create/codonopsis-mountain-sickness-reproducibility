options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
module_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
suppressPackageStartupMessages({library(glmnet); library(randomForest); library(e1071)})
source(file.path(script_dir, "final_ml_utils.R"))

out <- Sys.getenv("ML_OUTPUT_DIR", file.path(module_dir, "reproduced_results"))
observed <- readRDS(file.path(out, "final_ML_observed_fit.rds"))
old_perm <- readRDS(file.path(module_dir, "data", "dependencies",
                              "revised_ML_permutation_results.rds"))
x <- observed$input$x; y <- observed$input$y; subjects <- observed$input$subjects
codes <- old_perm$GSE103940_codes
stopifnot(length(codes) == 499L, !anyDuplicated(codes), all(codes > 0L))
unique_subjects <- sort(unique(subjects))

run_one <- function(j) {
  bits <- as.logical(as.integer(intToBits(codes[j]))[seq_along(unique_subjects)])
  yp <- paired_permute_labels(y, subjects, bits)
  fit <- run_final_loso(x, yp, subjects, 103940L, FALSE)
  perf <- final_performance(fit$predictions)
  perf$permutation_id <- sprintf("paired_perm_%04d", j)
  perf$swap_code <- codes[j]
  perf
}
cores <- max(1L, as.integer(Sys.getenv("ML_FINAL_PERMUTATION_CORES", "16")))
message("Running 499 paired assignments with ", cores, " processes")
if (cores == 1L) {
  null <- lapply(seq_along(codes), run_one)
} else {
  null <- parallel::mclapply(seq_along(codes), run_one, mc.cores = cores,
                             mc.preschedule = FALSE, mc.set.seed = FALSE)
}
null <- do.call(rbind, null)
write.csv(null, file.path(out, "GSE103940", "GSE103940_paired_permutation_null_metrics.csv"), row.names = FALSE)

summary_rows <- list(); k <- 0L
for (method in names(final_methods)) {
  obs <- observed$performance[observed$performance$method == method, ]
  z <- null[null$method == method, ]
  for (metric in c("ROC_AUC", "balanced_accuracy")) {
    k <- k + 1L
    value <- obs[[metric]]
    summary_rows[[k]] <- data.frame(
      method = method, metric = metric, observed = value,
      null_mean = mean(z[[metric]], na.rm = TRUE),
      null_median = median(z[[metric]], na.rm = TRUE),
      null_q95 = quantile(z[[metric]], .95, na.rm = TRUE),
      empirical_permutation_P = (1 + sum(z[[metric]] >= value, na.rm = TRUE)) /
        (nrow(z) + 1), permutations = nrow(z),
      permutation_scheme = "within-subject plain/high label swap")
  }
}
summary <- do.call(rbind, summary_rows)
write.csv(summary, file.path(out, "GSE103940", "GSE103940_paired_permutation_summary.csv"), row.names = FALSE)
saveRDS(list(null = null, summary = summary, codes = codes),
        file.path(out, "final_ML_paired_permutation_results.rds"))
capture.output(sessionInfo(), file = file.path(out, "final_ML_permutation_sessionInfo.txt"))
print(summary[summary$metric == "ROC_AUC", ])
