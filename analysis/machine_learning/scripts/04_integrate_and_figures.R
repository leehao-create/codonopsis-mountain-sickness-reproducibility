options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
module_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
suppressPackageStartupMessages(library(ggplot2))
out <- Sys.getenv("ML_OUTPUT_DIR", file.path(module_dir, "reproduced_results"))
figdir <- file.path(out, "fig5_components")
supdir <- file.path(out, "supplementary_GSE75665")
obs <- readRDS(file.path(out, "final_ML_observed_fit.rds"))
perm <- readRDS(file.path(out, "final_ML_paired_permutation_results.rds"))
stab <- read.csv(file.path(out, "ML_final_feature_stability.csv"), check.names = FALSE)
paired <- obs$paired
methods <- c("LASSO", "ElasticNet", "RF", "SVM", "Ensemble")
cols <- c(LASSO = "#0072B2", ElasticNet = "#009E73", RF = "#D55E00",
          SVM = "#CC79A7", Ensemble = "#222222")

save_plot <- function(plot, stem, width, height) {
  ggsave(file.path(figdir, paste0(stem, ".pdf")), plot, width = width,
         height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(figdir, paste0(stem, ".png")), plot, width = width,
         height = height, units = "in", dpi = 300)
}

# A: analysis workflow, deliberately schematic rather than a result panel.
workflow <- data.frame(
  x = 1:6, y = 1,
  label = c("54 predefined\ncandidate targets",
            "39 expression-eligible\npredictors",
            "Outer LOSO\n11 subjects",
            "Training-fold only\npreprocessing + tuning",
            "OOF prediction +\npaired permutation",
            "Stability + subject\nbootstrap uncertainty")
)
pA <- ggplot(workflow, aes(x, y)) +
  geom_segment(data = workflow[-nrow(workflow), ],
               aes(x = x + .33, xend = x + .67, yend = y),
               arrow = arrow(length = unit(0.12, "inches")), linewidth = .5,
               colour = "#555555") +
  geom_label(aes(label = label), size = 3.2, linewidth = .35,
             label.padding = unit(.35, "lines"), fill = "white") +
  annotate("text", x = 3.5, y = .67,
           label = "LASSO | Elastic Net | RF | linear SVM-RFE | equal-weight ensemble",
           size = 3.1, colour = "#333333") +
  coord_cartesian(xlim = c(.5, 6.5), ylim = c(.55, 1.25), clip = "off") +
  theme_void() + theme(plot.margin = margin(10, 20, 10, 20))
save_plot(pA, "Fig5A_final_ML_workflow", 12, 2.7)
write.csv(workflow, file.path(figdir, "Fig5A_workflow_data.csv"), row.names = FALSE)

roc_points <- function(y, score, method) {
  thresholds <- c(Inf, sort(unique(score), decreasing = TRUE), -Inf)
  do.call(rbind, lapply(thresholds, function(th) {
    pred <- score >= th
    data.frame(method = method,
               FPR = mean(pred[y == 0L]), TPR = mean(pred[y == 1L]))
  }))
}
pred <- obs$fit$predictions
roc <- do.call(rbind, lapply(methods, function(m) {
  column <- c(LASSO = "LASSO_probability", ElasticNet = "ElasticNet_probability",
              RF = "RF_probability", SVM = "SVM_probability",
              Ensemble = "Ensemble_probability")[[m]]
  roc_points(pred$observed, pred[[column]], m)
}))
auc_labels <- setNames(sprintf("%s (AUC %.3f)", obs$performance$method,
                               obs$performance$ROC_AUC), obs$performance$method)
roc$method_label <- factor(auc_labels[roc$method], levels = auc_labels[methods])
pB <- ggplot(roc, aes(FPR, TPR, colour = method_label)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#999999") +
  geom_step(linewidth = .9) + coord_equal() +
  scale_colour_manual(values = unname(cols[methods]), name = NULL) +
  labs(x = "False-positive rate", y = "True-positive rate",
       title = "GSE103940 outer-LOSO out-of-fold ROC") +
  theme_classic(base_size = 11) + theme(legend.position = "right")
save_plot(pB, "Fig5B_GSE103940_OOF_ROC", 6.7, 5.0)
write.csv(roc, file.path(figdir, "Fig5B_OOF_ROC_coordinates.csv"), row.names = FALSE)

null_auc <- perm$null
obs_auc <- obs$performance[, c("method", "ROC_AUC")]
pC <- ggplot(null_auc, aes(ROC_AUC)) +
  geom_histogram(bins = 25, fill = "#B8C4CE", colour = "white") +
  geom_vline(data = obs_auc, aes(xintercept = ROC_AUC), colour = "#B2182B",
             linewidth = .8) +
  facet_wrap(~factor(method, levels = methods), ncol = 3) +
  labs(x = "AUC under paired-label permutation", y = "Count",
       title = "Observed AUC (red) versus structure-preserving null") +
  theme_classic(base_size = 10)
save_plot(pC, "Fig5C_GSE103940_permutation_distributions", 8, 5.5)

long_stab <- reshape(stab, varying = list(c("LASSO_frequency", "ElasticNet_frequency",
                                            "RF_frequency", "SVM_frequency")),
                     v.names = "frequency", timevar = "method",
                     times = c("LASSO", "Elastic Net", "RF", "SVM"),
                     direction = "long")
long_stab <- long_stab[!is.na(long_stab$frequency), ]
long_stab$method <- factor(long_stab$method,
                           levels = c("LASSO", "Elastic Net", "RF", "SVM"))
gene_order <- aggregate(frequency ~ gene, long_stab, mean)
gene_order <- gene_order$gene[order(gene_order$frequency)]
long_stab$gene <- factor(long_stab$gene, levels = gene_order)
pD <- ggplot(long_stab, aes(method, gene, fill = frequency)) +
  geom_tile(colour = "white", linewidth = .25) +
  geom_text(data = long_stab[long_stab$frequency >= .60, ], label = "*", size = 3) +
  scale_fill_gradientn(colours = c("#F7FBFF", "#6BAED6", "#08306B"),
                       limits = c(0, 1), name = "Selection\nfrequency") +
  labs(x = NULL, y = NULL, title = "Candidate feature-selection stability",
       subtitle = "* predefined stable support: frequency >= 0.60") +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(pD, "Fig5D_GSE103940_feature_stability_heatmap", 6.5, 8.5)
write.csv(long_stab, file.path(figdir, "Fig5D_feature_stability_long.csv"), row.names = FALSE)

consensus <- stab[stab$number_of_methods_supported > 0L, ]
consensus$mean_frequency <- rowMeans(consensus[, 2:5], na.rm = TRUE)
consensus <- consensus[order(consensus$number_of_methods_supported,
                             consensus$mean_frequency), ]
consensus$gene <- factor(consensus$gene, levels = consensus$gene)
pE <- ggplot(consensus, aes(number_of_methods_supported, gene,
                            size = mean_frequency, colour = number_of_methods_supported)) +
  geom_point() +
  scale_colour_gradient(low = "#56B4E9", high = "#B2182B", limits = c(1, 4)) +
  scale_x_continuous(breaks = 1:4, limits = c(.7, 4.3)) +
  labs(x = "Methods with selection frequency >= 0.60", y = NULL,
       size = "Mean frequency", colour = "Methods",
       title = "Consensus candidate prioritization") +
  theme_classic(base_size = 10)
save_plot(pE, "Fig5E_GSE103940_consensus_candidates", 6.8, 5.8)
write.csv(consensus, file.path(figdir, "Fig5E_consensus_candidate_data.csv"), row.names = FALSE)

paired$method <- factor(paired$method, levels = methods)
pF <- ggplot(paired, aes(subject_id, score_difference_high_minus_plain,
                         colour = method, group = method)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "#777777") +
  geom_line(linewidth = .45, alpha = .55) + geom_point(size = 1.8) +
  scale_colour_manual(values = cols, name = NULL) +
  labs(x = "Held-out subject", y = "Predicted P(high): high-altitude - plain",
       title = "Paired out-of-fold score differences") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
save_plot(pF, "Fig5F_GSE103940_paired_score_differences", 8.2, 4.8)

# Frozen GSE75665 exploratory files are distributed under results/
# supplementary_GSE75665 and are intentionally not regenerated here.

capture.output(sessionInfo(), file = file.path(out, "final_ML_figure_sessionInfo.txt"))
