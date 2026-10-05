options(stringsAsFactors = FALSE)
set.seed(75665)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
module_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
suppressPackageStartupMessages({
  library(WGCNA)
  library(limma)
  library(lmerTest)
  library(ggplot2)
})
allowWGCNAThreads(nThreads = 4)

input_file <- file.path(module_dir, "data", "processed", "GSE75665_WGCNA_input_FPKM.csv")
metadata_file <- file.path(module_dir, "data", "metadata", "GSE75665_subject_pairing.csv")
target_file <- file.path(module_dir, "data", "reference", "Codonopsis_54_overlap_targets.csv")
out <- Sys.getenv("WGCNA_OUTPUT_DIR", file.path(module_dir, "reproduced_results"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(input_file), file.exists(metadata_file), file.exists(target_file))

expr_df <- read.csv(input_file, row.names = 1, check.names = FALSE)
expr_df[] <- lapply(expr_df, as.numeric)
expr <- as.matrix(expr_df)
meta <- read.csv(metadata_file, check.names = FALSE)
targets <- unique(read.csv(target_file, check.names = FALSE)[, 1])

stopifnot(!anyDuplicated(rownames(expr)), !anyDuplicated(colnames(expr)),
          !anyDuplicated(meta$gsm), identical(colnames(expr), meta$gsm),
          all(is.finite(expr)), all(expr >= 0), all(table(meta$subject_id) == 2L),
          length(targets) == 54L)

meta$AMS_status <- factor(meta$ams_status, levels = c("non_AMS", "AMS"))
meta$altitude_status <- factor(meta$altitude_status, levels = c("plain", "high_altitude"))
meta$subject_id <- factor(meta$subject_id)

log_expr <- log2(expr + 1)
valid_gene_symbol <- nzchar(rownames(expr)) & !is.na(rownames(expr))
eligible <- rowSums(expr >= 1) >= 5 & valid_gene_symbol
mad_log <- apply(log_expr, 1, mad, na.rm = TRUE)
nonzero_eligible <- eligible & mad_log > 0
q50 <- unname(quantile(mad_log[nonzero_eligible], 0.50))
q75 <- unname(quantile(mad_log[nonzero_eligible], 0.75))
filter_sets <- list(
  primary_MAD_upper50 = rownames(expr)[nonzero_eligible & mad_log >= q50],
  sensitivity_MAD_upper25 = rownames(expr)[nonzero_eligible & mad_log >= q75],
  sensitivity_all_expression_eligible = rownames(expr)[eligible]
)

filter_summary <- data.frame(
  stage = c("input", "FPKM >= 1 in >= 5/20 samples", "nonzero log2-scale MAD",
            "primary: upper 50% MAD", "sensitivity: upper 25% MAD",
            "sensitivity: all expression-eligible"),
  genes = c(nrow(expr), sum(eligible), sum(nonzero_eligible),
            length(filter_sets$primary_MAD_upper50),
            length(filter_sets$sensitivity_MAD_upper25),
            length(filter_sets$sensitivity_all_expression_eligible)),
  rule = c("none", "phenotype-blind expression eligibility",
           "remove zero-variance genes", paste0("MAD >= median = ", signif(q50, 6)),
           paste0("MAD >= 75th percentile = ", signif(q75, 6)),
           "no variance quantile filter")
)
write.csv(filter_summary, file.path(out, "WGCNA_gene_filter_summary.csv"), row.names = FALSE)

distribution_summary <- data.frame(
  metric = c("FPKM minimum", "FPKM Q1", "FPKM median", "FPKM Q3", "FPKM maximum",
             "log2(FPKM+1) minimum", "log2(FPKM+1) Q1", "log2(FPKM+1) median",
             "log2(FPKM+1) Q3", "log2(FPKM+1) maximum", "zero fraction"),
  value = c(min(expr), quantile(expr, 0.25), median(expr), quantile(expr, 0.75), max(expr),
            min(log_expr), quantile(log_expr, 0.25), median(log_expr),
            quantile(log_expr, 0.75), max(log_expr), mean(expr == 0))
)
write.csv(distribution_summary, file.path(out, "WGCNA_expression_distribution_summary.csv"),
          row.names = FALSE)

mad_plot_data <- data.frame(MAD = mad_log[eligible])
mad_plot <- ggplot(mad_plot_data, aes(MAD)) +
  geom_histogram(bins = 70, fill = "#4C78A8", color = "white", linewidth = 0.15) +
  geom_vline(xintercept = c(q50, q75), color = c("#D95F02", "#1B9E77"),
             linetype = c(1, 2)) +
  coord_cartesian(xlim = quantile(mad_plot_data$MAD, c(0, 0.99))) +
  labs(x = "MAD of log2(FPKM + 1)", y = "Eligible genes",
       title = "GSE75665 phenotype-blind variance filtering") +
  theme_classic(base_size = 12)
ggsave(file.path(out, "WGCNA_MAD_filter_distribution.pdf"), mad_plot, width = 7, height = 5)
ggsave(file.path(out, "WGCNA_MAD_filter_distribution.png"), mad_plot,
       width = 7, height = 5, dpi = 300)

primary_genes <- filter_sets$primary_MAD_upper50
datExpr <- t(log_expr[primary_genes, , drop = FALSE])
rownames(datExpr) <- meta$gsm
gsg <- goodSamplesGenes(datExpr, verbose = 3)
stopifnot(gsg$allOK)

# Prespecified objective sample-outlier rule.
sample_cor <- stats::cor(t(datExpr), use = "pairwise.complete.obs")
sample_adj <- (1 + sample_cor) / 2
diag(sample_adj) <- 0
sample_k <- rowSums(sample_adj)
sample_Zk <- as.numeric(scale(sample_k))
sample_dist <- as.matrix(dist(datExpr))
mean_dist <- rowSums(sample_dist) / (nrow(sample_dist) - 1)
robust_dist_Z <- (mean_dist - median(mean_dist)) / mad(mean_dist)
outlier_flag <- sample_Zk < -2.5 & robust_dist_Z > 3.5
sample_qc <- data.frame(
  sample_id = meta$gsm, subject_id = meta$subject_id, AMS_status = meta$AMS_status,
  altitude_status = meta$altitude_status, sample_connectivity = sample_k,
  connectivity_Z = sample_Zk, mean_inter_sample_distance = mean_dist,
  robust_distance_Z = robust_dist_Z, outlier = outlier_flag,
  decision = ifelse(outlier_flag, "remove by prespecified joint rule", "retain")
)
write.csv(sample_qc, file.path(out, "WGCNA_sample_QC.csv"), row.names = FALSE)

sample_tree <- hclust(dist(datExpr), method = "average")
pdf(file.path(out, "WGCNA_sample_clustering.pdf"), width = 12, height = 7)
plot(sample_tree,
     labels = paste(meta$gsm, meta$AMS_status, meta$altitude_status, sep = " | "),
     main = "GSE75665 sample clustering", sub = "", xlab = "", cex = 0.7)
dev.off()

if (any(outlier_flag)) {
  datExpr <- datExpr[!outlier_flag, , drop = FALSE]
  meta <- meta[!outlier_flag, , drop = FALSE]
  writeLines(paste(sum(outlier_flag), "sample(s) removed by the prespecified joint rule."),
             file.path(out, "WGCNA_sample_outlier_decision.txt"))
} else {
  writeLines("No obvious sample outliers were identified.",
             file.path(out, "WGCNA_sample_outlier_decision.txt"))
}
stopifnot(identical(rownames(datExpr), meta$gsm))

powers <- 1:30
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed",
                         corFnc = "cor", corOptions = list(use = "p"),
                         verbose = 5, blockSize = ncol(datExpr))
fit <- as.data.frame(sft$fitIndices)
required_cols <- c("Power", "SFT.R.sq", "slope", "mean.k.")
stopifnot(all(required_cols %in% colnames(fit)))
fit$negative_slope <- fit$slope < 0
fit$meets_R2_0_85 <- fit$SFT.R.sq >= 0.85
fit$nondegenerate_connectivity <- fit$mean.k. >= 1
fit$meets_prespecified_rule <- fit$negative_slope & fit$meets_R2_0_85 &
  fit$nondegenerate_connectivity
fit$within_prespecified_range <- TRUE
qualifying <- fit$Power[fit$meets_prespecified_rule]
if (length(qualifying)) {
  selected_power <- min(qualifying)
  power_reason <- paste0("Lowest power with R2 >= 0.85, negative slope, and mean connectivity >= 1: beta=",
                         selected_power)
  threshold_met <- TRUE
} else {
  compromise <- fit[fit$negative_slope & fit$nondegenerate_connectivity, ]
  stopifnot(nrow(compromise) > 0L)
  selected_power <- compromise$Power[which.max(compromise$SFT.R.sq)]
  power_reason <- paste0("No power met R2 >= 0.85; selected maximum R2 among negative-slope powers with mean connectivity >= 1: beta=",
                         selected_power)
  threshold_met <- FALSE
}
fit$selected <- fit$Power == selected_power
fit_out <- fit
if (!threshold_met) {
  sft_extension <- pickSoftThreshold(datExpr, powerVector = 31:50, networkType = "signed",
                                     corFnc = "cor", corOptions = list(use = "p"),
                                     verbose = 0, blockSize = ncol(datExpr))
  extension <- as.data.frame(sft_extension$fitIndices)
  extension$negative_slope <- extension$slope < 0
  extension$meets_R2_0_85 <- extension$SFT.R.sq >= 0.85
  extension$nondegenerate_connectivity <- extension$mean.k. >= 1
  extension$meets_prespecified_rule <- extension$negative_slope &
    extension$meets_R2_0_85 & extension$nondegenerate_connectivity
  extension$within_prespecified_range <- FALSE
  extension$selected <- FALSE
  fit_out <- rbind(fit, extension[, colnames(fit)])
  power_reason <- paste0(power_reason, ". Diagnostic extension through beta=50 also failed R2=0.85 (maximum extended R2=",
                         signif(max(extension$SFT.R.sq), 4), ")")
}
write.csv(fit_out, file.path(out, "WGCNA_soft_threshold_candidates.csv"), row.names = FALSE)
write.csv(data.frame(selected_power = selected_power, R2_threshold_met = threshold_met,
                     decision_rule = power_reason),
          file.path(out, "WGCNA_soft_threshold_decision.csv"), row.names = FALSE)

pdf(file.path(out, "WGCNA_soft_threshold.pdf"), width = 11, height = 5)
par(mfrow = c(1, 2))
plot(fit_out$Power, fit_out$SFT.R.sq, type = "n", xlab = "Soft-threshold power",
     ylab = "Scale-free topology fit R2", main = "Scale independence")
text(fit_out$Power, fit_out$SFT.R.sq, labels = fit_out$Power,
     col = ifelse(fit_out$selected, "red3",
                  ifelse(fit_out$within_prespecified_range, "grey30", "grey65")), cex = 0.65)
abline(h = 0.85, col = "red3", lty = 2)
plot(fit_out$Power, fit_out$mean.k., type = "n", xlab = "Soft-threshold power",
     ylab = "Mean connectivity", main = "Mean connectivity")
text(fit_out$Power, fit_out$mean.k., labels = fit_out$Power,
     col = ifelse(fit_out$selected, "red3",
                  ifelse(fit_out$within_prespecified_range, "grey30", "grey65")), cex = 0.65)
abline(h = 1, col = "grey50", lty = 2)
dev.off()

network_args <- list(
  power = selected_power, networkType = "signed", TOMType = "signed",
  corType = "pearson", minModuleSize = 20, deepSplit = 2,
  mergeCutHeight = 0.25, reassignThreshold = 0,
  pamRespectsDendro = FALSE, numericLabels = FALSE,
  maxBlockSize = 6000, saveTOMs = FALSE, verbose = 3
)
net <- do.call(blockwiseModules, c(list(datExpr = datExpr), network_args))
module_colors <- net$colors
names(module_colors) <- colnames(datExpr)
MEs <- orderMEs(net$MEs)
# Grey denotes unassigned genes, not a biological module.
MEs <- MEs[, colnames(MEs) != "MEgrey", drop = FALSE]
stopifnot(nrow(MEs) == nrow(meta))

module_assignment <- data.frame(gene = names(module_colors), module = unname(module_colors))
write.csv(module_assignment, file.path(out, "WGCNA_module_assignments.csv"), row.names = FALSE)
write.csv(data.frame(sample_id = meta$gsm, MEs, check.names = FALSE),
          file.path(out, "WGCNA_module_eigengenes.csv"), row.names = FALSE)
write.csv(as.data.frame(table(module_colors)), file.path(out, "WGCNA_module_sizes.csv"),
          row.names = FALSE)

pdf(file.path(out, "WGCNA_gene_dendrogram_modules.pdf"), width = 13, height = 7)
for (block in seq_along(net$dendrograms)) {
  block_genes <- net$blockGenes[[block]]
  plotDendroAndColors(net$dendrograms[[block]], module_colors[block_genes],
                      "Signed modules", dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = paste("Gene dendrogram and modules - block", block))
}
dev.off()

# Subject-aware mixed models for every module eigengene.
fixed_terms <- c("AMS_statusAMS", "altitude_statushigh_altitude",
                 "AMS_statusAMS:altitude_statushigh_altitude")
mixed_rows <- list()
for (module in colnames(MEs)) {
  model_df <- data.frame(ME = MEs[, module], subject_id = meta$subject_id,
                         AMS_status = meta$AMS_status,
                         altitude_status = meta$altitude_status)
  model <- lmerTest::lmer(ME ~ AMS_status * altitude_status + (1 | subject_id),
                          data = model_df, REML = TRUE)
  co <- coef(summary(model))
  stopifnot(all(fixed_terms %in% rownames(co)))
  singular <- lme4::isSingular(model, tol = 1e-5)
  for (term in fixed_terms) {
    mixed_rows[[length(mixed_rows) + 1L]] <- data.frame(
      module = module, effect = term, coefficient = co[term, "Estimate"],
      standard_error = co[term, "Std. Error"], df = co[term, "df"],
      t_value = co[term, "t value"], P_value = co[term, "Pr(>|t|)"],
      standardized_effect = co[term, "Estimate"] / sd(model_df$ME),
      singular_fit = singular
    )
  }
}
mixed <- do.call(rbind, mixed_rows)
mixed$BH_FDR_global <- p.adjust(mixed$P_value, method = "BH")
mixed$BH_FDR_within_effect <- ave(mixed$P_value, mixed$effect,
                                  FUN = function(z) p.adjust(z, method = "BH"))
mixed <- mixed[order(mixed$P_value), ]
write.csv(mixed, file.path(out, "WGCNA_module_trait_mixed_effects.csv"), row.names = FALSE)

interaction_term <- "AMS_statusAMS:altitude_statushigh_altitude"
sig_modules <- sub("^ME", "", mixed$module[mixed$effect == interaction_term &
                                             mixed$BH_FDR_global < 0.05])

# Equivalent subject-level change-score sensitivity for the interaction.
delta_ME <- do.call(rbind, lapply(levels(meta$subject_id), function(s) {
  idx <- which(meta$subject_id == s)
  stopifnot(length(idx) == 2L)
  high <- idx[meta$altitude_status[idx] == "high_altitude"]
  plain <- idx[meta$altitude_status[idx] == "plain"]
  stopifnot(length(high) == 1L, length(plain) == 1L)
  MEs[high, , drop = FALSE] - MEs[plain, , drop = FALSE]
}))
rownames(delta_ME) <- levels(meta$subject_id)
subject_status <- meta$AMS_status[match(rownames(delta_ME), meta$subject_id)]
change_rows <- lapply(colnames(delta_ME), function(module) {
  x <- delta_ME[, module]
  a <- x[subject_status == "AMS"]
  n <- x[subject_status == "non_AMS"]
  test <- t.test(a, n)
  pooled_sd <- sqrt(((length(a) - 1) * var(a) + (length(n) - 1) * var(n)) /
                      (length(a) + length(n) - 2))
  data.frame(module = module, mean_change_AMS = mean(a),
             mean_change_non_AMS = mean(n), difference_in_change = mean(a) - mean(n),
             Cohen_d = (mean(a) - mean(n)) / pooled_sd, P_value = test$p.value)
})
change_score <- do.call(rbind, change_rows)
change_score$BH_FDR <- p.adjust(change_score$P_value, method = "BH")
write.csv(change_score[order(change_score$P_value), ],
          file.path(out, "WGCNA_module_interaction_change_score.csv"), row.names = FALSE)

# Module-trait result figure: coefficients with global FDR annotations.
mixed$module_label <- sub("^ME", "", mixed$module)
mixed$effect_label <- factor(mixed$effect, levels = fixed_terms,
                             labels = c("AMS at plain", "Altitude in non-AMS", "AMS x altitude"))
mixed$label <- paste0(sprintf("%.2f", mixed$coefficient), "\nFDR=",
                      formatC(mixed$BH_FDR_global, format = "g", digits = 2))
trait_plot <- ggplot(mixed, aes(effect_label, module_label, fill = coefficient)) +
  geom_tile(color = "white") + geom_text(aes(label = label), size = 2.5) +
  scale_fill_gradient2(low = "#2878B5", mid = "white", high = "#C33C3C") +
  labs(x = NULL, y = "Module", fill = "Coefficient",
       title = "Subject-aware module-phenotype models") +
  theme_minimal(base_size = 10) + theme(panel.grid = element_blank())
ggsave(file.path(out, "WGCNA_module_trait_mixed_effects.pdf"), trait_plot,
       width = 9, height = max(5, 0.38 * length(unique(mixed$module_label))))
ggsave(file.path(out, "WGCNA_module_trait_mixed_effects.png"), trait_plot,
       width = 9, height = max(5, 0.38 * length(unique(mixed$module_label))), dpi = 300)

# Formal enrichment of the 54 prior targets in every module.
universe <- names(module_colors)
target_universe <- intersect(targets, universe)
enrichment <- do.call(rbind, lapply(sort(setdiff(unique(module_colors), "grey")), function(module) {
  genes <- names(module_colors)[module_colors == module]
  k <- sum(genes %in% target_universe)
  n <- length(genes)
  K <- length(target_universe)
  N <- length(universe)
  ft <- fisher.test(matrix(c(k, n - k, K - k, N - n - K + k),
                           nrow = 2, byrow = TRUE), alternative = "greater")
  data.frame(module = module, module_size = n, overlapping_targets = k,
             targets_in_network_universe = K, odds_ratio = unname(ft$estimate),
             P_value = ft$p.value,
             target_genes = paste(sort(intersect(genes, target_universe)), collapse = ";"))
}))
enrichment$BH_FDR <- p.adjust(enrichment$P_value, method = "BH")
enrichment <- enrichment[order(enrichment$P_value), ]
write.csv(enrichment, file.path(out, "WGCNA_Codonopsis_target_enrichment.csv"),
          row.names = FALSE)

# Subject-aware gene significance from per-subject altitude change.
delta_gene <- do.call(rbind, lapply(levels(meta$subject_id), function(s) {
  idx <- which(meta$subject_id == s)
  high <- idx[meta$altitude_status[idx] == "high_altitude"]
  plain <- idx[meta$altitude_status[idx] == "plain"]
  datExpr[high, , drop = FALSE] - datExpr[plain, , drop = FALSE]
}))
rownames(delta_gene) <- levels(meta$subject_id)
gene_design <- model.matrix(~ subject_status)
fit_gene <- eBayes(lmFit(t(delta_gene), gene_design), robust = TRUE)
gene_gs <- topTable(fit_gene, coef = "subject_statusAMS", number = Inf,
                    adjust.method = "BH", sort.by = "none")
gene_gs$gene <- rownames(gene_gs)

# Module membership and intramodular connectivity.
kME_matrix <- stats::cor(datExpr, MEs, use = "pairwise.complete.obs")
kME_assigned <- vapply(names(module_colors), function(g) {
  module_col <- paste0("ME", module_colors[[g]])
  if (module_col %in% colnames(kME_matrix)) kME_matrix[g, module_col] else NA_real_
}, numeric(1))
adj <- adjacency(datExpr, power = selected_power, type = "signed", corFnc = "cor",
                 corOptions = list(use = "p"))
connectivity <- intramodularConnectivity(adj, module_colors, scaleByMax = FALSE)
rm(adj)
gc()

metrics <- data.frame(
  gene = names(module_colors), module = unname(module_colors),
  kME = unname(kME_assigned),
  kTotal = connectivity[names(module_colors), "kTotal"],
  kWithin = connectivity[names(module_colors), "kWithin"],
  kOut = connectivity[names(module_colors), "kOut"],
  stringsAsFactors = FALSE
)
metrics$kWithin_90th_percentile <- ave(metrics$kWithin, metrics$module,
                                      FUN = function(z) quantile(z, 0.90, na.rm = TRUE))
metrics$gene_interaction_logFC <- gene_gs[metrics$gene, "logFC"]
metrics$gene_interaction_t <- gene_gs[metrics$gene, "t"]
metrics$gene_interaction_P <- gene_gs[metrics$gene, "P.Value"]
metrics$gene_interaction_FDR <- gene_gs[metrics$gene, "adj.P.Val"]
metrics$GS_abs_t <- abs(metrics$gene_interaction_t)
metrics$module_interaction_FDR_significant <- metrics$module %in% sig_modules
metrics$hub_abs_kME_ge_0_80 <- abs(metrics$kME) >= 0.80
metrics$hub_kWithin_top10pct <- metrics$kWithin >= metrics$kWithin_90th_percentile
metrics$hub_gene_interaction_FDR_lt_0_05 <- metrics$gene_interaction_FDR < 0.05
metrics$hub_gene <- metrics$module_interaction_FDR_significant &
  metrics$hub_abs_kME_ge_0_80 & metrics$hub_kWithin_top10pct &
  metrics$hub_gene_interaction_FDR_lt_0_05
write.csv(metrics, file.path(out, "WGCNA_gene_module_metrics.csv"), row.names = FALSE)
write.csv(metrics[metrics$hub_gene, ], file.path(out, "WGCNA_hub_genes.csv"), row.names = FALSE)

candidates <- c("NR3C2", "MMP9", "SNCA", "HBB", "ABCB1", "BCL2", "CA2")
candidate_metrics <- merge(data.frame(gene = candidates), metrics, by = "gene", all.x = TRUE,
                           sort = FALSE)
candidate_metrics$in_primary_network <- !is.na(candidate_metrics$module)
write.csv(candidate_metrics, file.path(out, "WGCNA_candidate_gene_metrics.csv"), row.names = FALSE)
mmp9 <- candidate_metrics[candidate_metrics$gene == "MMP9", ]
if (nrow(mmp9) == 1L && isTRUE(mmp9$hub_gene)) {
  mmp9_statement <- "MMP9 satisfies the prespecified revised WGCNA hub criteria."
} else {
  mmp9_statement <- "MMP9 is not supported as a WGCNA hub gene under the revised analysis."
}
writeLines(mmp9_statement, file.path(out, "MMP9_WGCNA_conclusion.txt"))

# Network stability across neighboring powers and prespecified filter variants.
adjusted_rand <- function(x, y) {
  tab <- table(x, y)
  n <- sum(tab)
  choose2 <- function(z) z * (z - 1) / 2
  a <- sum(choose2(tab))
  row_sum <- sum(choose2(rowSums(tab)))
  col_sum <- sum(choose2(colSums(tab)))
  total <- choose2(n)
  expected <- row_sum * col_sum / total
  denom <- 0.5 * (row_sum + col_sum) - expected
  if (denom == 0) return(NA_real_)
  (a - expected) / denom
}

run_sensitivity_network <- function(genes, power, name) {
  d <- t(log_expr[genes, meta$gsm, drop = FALSE])
  n <- blockwiseModules(d, power = power, networkType = "signed", TOMType = "signed",
                        corType = "pearson", minModuleSize = 20, deepSplit = 2,
                        mergeCutHeight = 0.25, reassignThreshold = 0,
                        pamRespectsDendro = FALSE, numericLabels = FALSE,
                        maxBlockSize = 5000, saveTOMs = FALSE, verbose = 1)
  colors <- n$colors
  names(colors) <- colnames(d)
  write.csv(data.frame(gene = names(colors), module = unname(colors)),
            file.path(out, paste0("WGCNA_sensitivity_assignments_", name, ".csv")),
            row.names = FALSE)
  common <- intersect(names(module_colors), names(colors))
  common_non_grey <- common[module_colors[common] != "grey" & colors[common] != "grey"]
  data.frame(
    scenario = name, power = power, genes = length(colors),
    modules_excluding_grey = length(setdiff(unique(colors), "grey")),
    grey_fraction = mean(colors == "grey"),
    largest_module_fraction = max(table(colors[colors != "grey"])) / length(colors),
    common_genes_with_primary = length(common),
    adjusted_Rand_all_common = adjusted_rand(module_colors[common], colors[common]),
    adjusted_Rand_non_grey = if (length(common_non_grey) > 1L)
      adjusted_rand(module_colors[common_non_grey], colors[common_non_grey]) else NA_real_
  )
}

primary_summary <- data.frame(
  scenario = "primary", power = selected_power, genes = length(module_colors),
  modules_excluding_grey = length(setdiff(unique(module_colors), "grey")),
  grey_fraction = mean(module_colors == "grey"),
  largest_module_fraction = max(table(module_colors[module_colors != "grey"])) /
    length(module_colors), common_genes_with_primary = length(module_colors),
  adjusted_Rand_all_common = 1, adjusted_Rand_non_grey = 1
)
sensitivity_rows <- list(primary_summary)
neighbor_powers <- unique(pmax(1, pmin(30, c(selected_power - 1, selected_power + 1, 16))))
neighbor_powers <- setdiff(neighbor_powers, selected_power)
for (p in neighbor_powers) {
  sensitivity_rows[[length(sensitivity_rows) + 1L]] <- run_sensitivity_network(
    primary_genes, p, paste0("primary_filter_beta_", p))
}
sensitivity_rows[[length(sensitivity_rows) + 1L]] <- run_sensitivity_network(
  filter_sets$sensitivity_MAD_upper25, selected_power, "MAD_upper25")
sensitivity_rows[[length(sensitivity_rows) + 1L]] <- run_sensitivity_network(
  filter_sets$sensitivity_all_expression_eligible, selected_power,
  "all_expression_eligible")
sensitivity <- do.call(rbind, sensitivity_rows)
write.csv(sensitivity, file.path(out, "WGCNA_network_sensitivity.csv"), row.names = FALSE)

config <- data.frame(
  item = c("input", "scale", "transformation", "expression_filter", "primary_MAD_filter",
           "sample_outlier_rule", "selected_power", "power_decision", "networkType", "TOMType",
           "correlation", "minModuleSize", "deepSplit", "mergeCutHeight", "hub_criteria", "seed"),
  value = c(input_file, "FPKM", "log2(FPKM + 1)", "FPKM >= 1 in >= 5/20 samples",
            "upper 50% nonzero MAD among expression-eligible genes",
            "connectivity Z < -2.5 AND robust distance Z > 3.5",
            selected_power, power_reason, "signed", "signed", "Pearson", "20", "2", "0.25",
            "interaction module global FDR<0.05; |kME|>=0.80; kWithin top10%; gene interaction FDR<0.05",
            "75665")
)
write.csv(config, file.path(out, "WGCNA_analysis_config.csv"), row.names = FALSE)
write.csv(data.frame(path = c(input_file, metadata_file, target_file),
                     md5 = unname(tools::md5sum(c(input_file, metadata_file, target_file)))),
          file.path(out, "WGCNA_input_md5.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "WGCNA_sessionInfo.txt"))
saveRDS(list(network = net, selected_power = selected_power, metadata = meta,
             filter_sets = filter_sets, mixed_effects = mixed),
        file.path(out, "WGCNA_primary_network.rds"))

summary_lines <- c(
  paste0("Selected beta: ", selected_power),
  paste0("Power decision: ", power_reason),
  paste0("Samples retained: ", nrow(meta), "/20"),
  paste0("Primary network genes: ", ncol(datExpr)),
  paste0("FDR-significant AMS x altitude modules: ",
         ifelse(length(sig_modules), paste(sig_modules, collapse = ", "), "none")),
  paste0("FDR-significant target-enriched modules: ",
         ifelse(any(enrichment$BH_FDR < 0.05),
                paste(enrichment$module[enrichment$BH_FDR < 0.05], collapse = ", "), "none")),
  mmp9_statement
)
writeLines(summary_lines, file.path(out, "WGCNA_key_results.txt"))
message(paste(summary_lines, collapse = "\n"))
