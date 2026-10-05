options(stringsAsFactors = FALSE)
set.seed(4103940)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
module_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(pheatmap)
  library(fgsea)
})

out <- Sys.getenv("DE_OUTPUT_DIR", file.path(module_dir, "reproduced_results", "GSE103940"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
count_file <- file.path(module_dir, "data", "processed", "GSE103940_raw_counts_GRCh38.p13_NCBI.tsv.gz")
meta_file <- file.path(module_dir, "data", "metadata", "GSE103940_subject_pairing.csv")
map_file <- Sys.getenv(
  "NCBI_GENE_MAP",
  file.path(module_dir, "data", "reference", "NCBI_GeneID_to_symbol.txt")
)
target_file <- file.path(module_dir, "..", "network_pharmacology", "data", "derived",
                         "Codonopsis_mountain_sickness_54_targets.csv")
phase3_file <- file.path(module_dir, "data", "dependencies", "GSE103940_NCBI_counts_voom_DEG_all.csv")
msigdb_dir <- Sys.getenv("MSIGDB_DIR", file.path(module_dir, "data", "gene_sets"))
gmt_files <- c(
  Hallmark = file.path(msigdb_dir, "h.all.v2025.1.Hs.symbols.gmt"),
  Reactome = file.path(msigdb_dir, "c2.cp.reactome.v2025.1.Hs.symbols.gmt"),
  GO_BP = file.path(msigdb_dir, "c5.go.bp.v2025.1.Hs.symbols.gmt")
)
if (!file.exists(map_file)) {
  stop("NCBI gene mapping is not redistributed. Run scripts/download_ncbi_gene_mapping.sh ",
       "or set NCBI_GENE_MAP; see the module README.")
}
stopifnot(dir.exists(out), all(file.exists(c(count_file, meta_file, target_file,
                                              phase3_file, gmt_files))))

save_plot <- function(plot, stem, width = 8, height = 6) {
  ggsave(file.path(out, paste0(stem, ".pdf")), plot, width = width, height = height)
  ggsave(file.path(out, paste0(stem, ".png")), plot, width = width, height = height,
         dpi = 300)
}

meta <- read.csv(meta_file, check.names = FALSE)
meta$condition <- factor(meta$condition, levels = c("plain", "high_altitude"))
meta$subject_id <- factor(meta$subject_id)
count_df <- read.delim(count_file, check.names = FALSE)
stopifnot(colnames(count_df)[1] == "GeneID", identical(colnames(count_df)[-1], meta$gsm),
          !anyDuplicated(count_df$GeneID), all(table(meta$subject_id) == 2L))
counts_id <- as.matrix(count_df[, -1, drop = FALSE])
storage.mode(counts_id) <- "integer"
rownames(counts_id) <- as.character(count_df$GeneID)
stopifnot(all(counts_id >= 0), all(is.finite(counts_id)))

gene_map <- read.delim(map_file, check.names = FALSE, quote = "", fill = TRUE)
colnames(gene_map)[1:2] <- c("GeneID", "gene")
gene_map$GeneID <- as.character(gene_map$GeneID)
gene_map$gene <- trimws(gene_map$gene)
gene_map <- gene_map[nzchar(gene_map$GeneID) & nzchar(gene_map$gene), c("GeneID", "gene")]
gene_map <- gene_map[!duplicated(gene_map$GeneID), ]
mapped_gene <- gene_map$gene[match(rownames(counts_id), gene_map$GeneID)]
mapped <- !is.na(mapped_gene) & nzchar(mapped_gene)
counts_symbol <- rowsum(counts_id[mapped, , drop = FALSE], group = mapped_gene[mapped],
                        reorder = FALSE)

design <- model.matrix(~ subject_id + condition, data = meta)
rownames(design) <- meta$gsm
stopifnot(qr(design)$rank == ncol(design), "conditionhigh_altitude" %in% colnames(design))
write.csv(data.frame(sample_id = rownames(design), design, check.names = FALSE),
          file.path(out, "GSE103940_final_design_matrix.csv"), row.names = FALSE)

y0 <- DGEList(counts = counts_symbol)
keep <- filterByExpr(y0, design = design)
y0 <- y0[keep, , keep.lib.sizes = FALSE]
globin_requested <- c("HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2", "HBM", "HBZ")
globin_in_matrix <- intersect(globin_requested, rownames(y0))
stopifnot(length(globin_in_matrix) >= 6L)

# Analysis A: standard TMM on all filtered genes.
y_a <- calcNormFactors(y0, method = "TMM")
v_a <- voom(y_a, design, plot = FALSE)
fit_a <- eBayes(lmFit(v_a, design), robust = TRUE)

# Analysis B: estimate TMM factors from non-globin genes while retaining original
# library sizes and then apply those factors to the complete filtered matrix.
y_b_factor <- calcNormFactors(y0[!rownames(y0) %in% globin_in_matrix, , keep.lib.sizes = TRUE],
                              method = "TMM")
y_b <- y0
y_b$samples$norm.factors <- y_b_factor$samples$norm.factors
v_b <- voom(y_b, design, plot = FALSE)
fit_b <- eBayes(lmFit(v_b, design), robust = TRUE)

make_table <- function(fit, voom_object, method) {
  tt <- topTable(fit, coef = "conditionhigh_altitude", number = Inf,
                 adjust.method = "BH", sort.by = "P")
  tt$gene <- rownames(tt)
  tt$method <- method
  tt$direction <- ifelse(tt$logFC > 0, "up_in_high_altitude",
                         ifelse(tt$logFC < 0, "down_in_high_altitude", "unchanged"))
  tt$passes_FDR_0_05 <- tt$adj.P.Val < 0.05
  tt$passes_primary_threshold <- tt$adj.P.Val < 0.05 & abs(tt$logFC) > 1.5
  tt <- tt[, c("gene", "method", setdiff(colnames(tt), c("gene", "method")))]
  tt
}

tt_a <- make_table(fit_a, v_a, "A_standard_TMM")
tt_b <- make_table(fit_b, v_b, "B_non_globin_factor_TMM")
write.csv(tt_a, file.path(out, "GSE103940_analysis_A_DEG_all.csv"), row.names = FALSE)
write.csv(tt_a[tt_a$passes_primary_threshold, ],
          file.path(out, "GSE103940_analysis_A_DEG_significant.csv"), row.names = FALSE)
write.csv(tt_b, file.path(out, "GSE103940_analysis_B_DEG_all.csv"), row.names = FALSE)
write.csv(tt_b[tt_b$passes_primary_threshold, ],
          file.path(out, "GSE103940_analysis_B_DEG_significant.csv"), row.names = FALSE)

# Analysis A must reproduce the frozen Phase 3 result exactly enough to exclude pipeline drift.
phase3 <- read.csv(phase3_file, check.names = FALSE)
phase3 <- phase3[match(tt_a$gene, phase3$gene), ]
stopifnot(identical(tt_a$gene, phase3$gene),
          max(abs(tt_a$logFC - phase3$logFC)) < 1e-10,
          max(abs(tt_a$t - phase3$t)) < 1e-10,
          max(abs(tt_a$adj.P.Val - phase3$adj.P.Val)) < 1e-10,
          sum(tt_a$passes_primary_threshold) == 296L)

comparison <- merge(
  tt_a[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val",
           "passes_FDR_0_05", "passes_primary_threshold")],
  tt_b[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val",
           "passes_FDR_0_05", "passes_primary_threshold")],
  by = "gene", suffixes = c("_A", "_B"), sort = FALSE)
comparison <- comparison[match(tt_a$gene, comparison$gene), ]
comparison$same_direction <- sign(comparison$logFC_A) == sign(comparison$logFC_B)
comparison$abs_logFC_difference <- abs(comparison$logFC_A - comparison$logFC_B)
comparison$rank_P_A <- rank(comparison$P.Value_A, ties.method = "average")
comparison$rank_P_B <- rank(comparison$P.Value_B, ties.method = "average")
write.csv(comparison, file.path(out, "GSE103940_A_vs_B_gene_comparison.csv"), row.names = FALSE)

jaccard <- function(a, b) {
  u <- union(a, b)
  if (!length(u)) return(1)
  length(intersect(a, b)) / length(u)
}
fdr_a <- tt_a$gene[tt_a$passes_FDR_0_05]
fdr_b <- tt_b$gene[tt_b$passes_FDR_0_05]
primary_a <- tt_a$gene[tt_a$passes_primary_threshold]
primary_b <- tt_b$gene[tt_b$passes_primary_threshold]
metrics <- data.frame(
  metric = c("tested_genes", "logFC_Pearson", "logFC_Spearman", "moderated_t_Pearson",
             "moderated_t_Spearman", "direction_concordance", "max_abs_logFC_difference",
             "FDR_A", "FDR_B", "FDR_overlap", "FDR_Jaccard",
             "primary_A", "primary_B", "primary_overlap", "primary_Jaccard"),
  value = c(nrow(comparison), cor(comparison$logFC_A, comparison$logFC_B),
            cor(comparison$logFC_A, comparison$logFC_B, method = "spearman"),
            cor(comparison$t_A, comparison$t_B),
            cor(comparison$t_A, comparison$t_B, method = "spearman"),
            mean(comparison$same_direction), max(comparison$abs_logFC_difference),
            length(fdr_a), length(fdr_b), length(intersect(fdr_a, fdr_b)), jaccard(fdr_a, fdr_b),
            length(primary_a), length(primary_b), length(intersect(primary_a, primary_b)),
            jaccard(primary_a, primary_b))
)
write.csv(metrics, file.path(out, "GSE103940_A_vs_B_consistency_metrics.csv"), row.names = FALSE)

counts_summary <- do.call(rbind, lapply(list(A_standard_TMM = tt_a,
                                             B_non_globin_factor_TMM = tt_b),
  function(tt) data.frame(
    tested = nrow(tt), FDR_only = sum(tt$passes_FDR_0_05),
    FDR_only_up = sum(tt$passes_FDR_0_05 & tt$logFC > 0),
    FDR_only_down = sum(tt$passes_FDR_0_05 & tt$logFC < 0),
    primary = sum(tt$passes_primary_threshold),
    primary_up = sum(tt$passes_primary_threshold & tt$logFC > 0),
    primary_down = sum(tt$passes_primary_threshold & tt$logFC < 0),
    median_logFC = median(tt$logFC))))
counts_summary$method <- rownames(counts_summary)
rownames(counts_summary) <- NULL
counts_summary <- counts_summary[, c("method", setdiff(colnames(counts_summary), "method"))]
write.csv(counts_summary, file.path(out, "GSE103940_A_vs_B_DEG_counts.csv"), row.names = FALSE)

norm_factors <- data.frame(
  sample_id = meta$gsm, subject_id = as.character(meta$subject_id),
  condition = as.character(meta$condition), raw_library_size = y0$samples$lib.size,
  norm_factor_A = y_a$samples$norm.factors,
  norm_factor_B = y_b$samples$norm.factors,
  effective_library_size_A = y_a$samples$lib.size * y_a$samples$norm.factors,
  effective_library_size_B = y_b$samples$lib.size * y_b$samples$norm.factors
)
norm_factors$factor_ratio_B_over_A <- norm_factors$norm_factor_B / norm_factors$norm_factor_A
write.csv(norm_factors, file.path(out, "GSE103940_normalization_factor_comparison.csv"),
          row.names = FALSE)

top_n <- 30L
top_a <- head(tt_a$gene, top_n)
top_b <- head(tt_b$gene, top_n)
top_union <- union(top_a, top_b)
top_table <- comparison[match(top_union, comparison$gene), ]
top_table$rank_A <- match(top_table$gene, tt_a$gene)
top_table$rank_B <- match(top_table$gene, tt_b$gene)
top_table$top30_A <- top_table$gene %in% top_a
top_table$top30_B <- top_table$gene %in% top_b
write.csv(top_table, file.path(out, "GSE103940_A_vs_B_top_genes.csv"), row.names = FALSE)

targets <- unique(trimws(read.csv(target_file, check.names = FALSE)[[1]]))
stopifnot(length(targets) == 54L)
target_results <- data.frame(gene = targets, expression_eligible = targets %in% tt_a$gene)
target_results <- merge(target_results, comparison, by = "gene", all.x = TRUE, sort = FALSE)
target_results <- target_results[match(targets, target_results$gene), ]
write.csv(target_results, file.path(out, "GSE103940_Codonopsis_targets_A_vs_B.csv"),
          row.names = FALSE)

# QC plots on the two normalized voom expression matrices.
voom_long <- rbind(
  data.frame(expression = as.vector(v_a$E), sample_id = rep(meta$gsm, each = nrow(v_a$E)),
             condition = rep(as.character(meta$condition), each = nrow(v_a$E)),
             method = "A: standard TMM"),
  data.frame(expression = as.vector(v_b$E), sample_id = rep(meta$gsm, each = nrow(v_b$E)),
             condition = rep(as.character(meta$condition), each = nrow(v_b$E)),
             method = "B: non-globin factor TMM")
)
density_plot <- ggplot(voom_long, aes(expression, group = sample_id, color = condition)) +
  geom_density(linewidth = 0.35, alpha = 0.5) + facet_wrap(~method, ncol = 1) +
  scale_color_manual(values = c(plain = "#2878B5", high_altitude = "#C33C3C")) +
  labs(x = "voom logCPM", y = "Density", color = NULL,
       title = "GSE103940 normalized expression density") +
  theme_classic(base_size = 10) + theme(legend.position = "top")
save_plot(density_plot, "GSE103940_A_vs_B_density", 9, 8)

box_plot <- ggplot(voom_long, aes(sample_id, expression, fill = condition)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.25) + facet_wrap(~method, ncol = 1) +
  scale_fill_manual(values = c(plain = "#6BAED6", high_altitude = "#E6550D")) +
  labs(x = NULL, y = "voom logCPM", fill = NULL,
       title = "GSE103940 normalized sample distributions") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "top")
save_plot(box_plot, "GSE103940_A_vs_B_boxplot", 11, 8)

make_pca <- function(e, method) {
  p <- prcomp(t(e), center = TRUE, scale. = FALSE)
  variance <- 100 * p$sdev^2 / sum(p$sdev^2)
  data.frame(sample_id = meta$gsm, subject_id = as.character(meta$subject_id),
             condition = as.character(meta$condition), PC1 = p$x[, 1], PC2 = p$x[, 2],
             PC1_percent = variance[1], PC2_percent = variance[2], method = method)
}
pca_df <- rbind(make_pca(v_a$E, "A: standard TMM"),
                make_pca(v_b$E, "B: non-globin factor TMM"))
write.csv(pca_df, file.path(out, "GSE103940_A_vs_B_PCA_scores.csv"), row.names = FALSE)
pca_plot <- ggplot(pca_df, aes(PC1, PC2, color = condition, group = subject_id)) +
  geom_line(color = "grey70", linewidth = 0.3) + geom_point(size = 2.2) +
  facet_wrap(~method, scales = "free") +
  scale_color_manual(values = c(plain = "#2878B5", high_altitude = "#C33C3C")) +
  labs(x = "PC1", y = "PC2", color = NULL,
       title = "GSE103940 PCA; lines connect paired samples") +
  theme_classic(base_size = 10) + theme(legend.position = "top")
save_plot(pca_plot, "GSE103940_A_vs_B_PCA", 11, 5.5)

ma_df <- rbind(tt_a, tt_b)
ma_df$method_label <- ifelse(ma_df$method == "A_standard_TMM", "A: standard TMM",
                             "B: non-globin factor TMM")
ma_df$category <- "not_primary"
ma_df$category[ma_df$passes_primary_threshold & ma_df$logFC > 0] <- "up"
ma_df$category[ma_df$passes_primary_threshold & ma_df$logFC < 0] <- "down"
ma_plot <- ggplot(ma_df, aes(AveExpr, logFC, color = category)) +
  geom_point(alpha = 0.5, size = 0.65) + facet_wrap(~method_label) +
  geom_hline(yintercept = 0, color = "grey30") +
  geom_hline(yintercept = c(-1.5, 1.5), linetype = 2, color = "grey50") +
  scale_color_manual(values = c(down = "#2878B5", not_primary = "#A7A9AC", up = "#C33C3C")) +
  labs(x = "Average voom logCPM", y = "log2FC (high - plain)", color = NULL,
       title = "GSE103940 MA plots") +
  theme_classic(base_size = 10) + theme(legend.position = "top")
save_plot(ma_plot, "GSE103940_A_vs_B_MA_plot", 11, 5.5)

fc_df <- rbind(data.frame(logFC = tt_a$logFC, method = "A: standard TMM"),
               data.frame(logFC = tt_b$logFC, method = "B: non-globin factor TMM"))
fc_plot <- ggplot(fc_df, aes(logFC, color = method)) +
  geom_density(linewidth = 0.8) + geom_vline(xintercept = 0, linetype = 2) +
  scale_color_manual(values = c("A: standard TMM" = "#2878B5",
                                "B: non-globin factor TMM" = "#C33C3C")) +
  labs(x = "log2FC (high altitude - plain)", y = "Density", color = NULL,
       title = "GSE103940 A/B logFC distributions") +
  theme_classic(base_size = 11) + theme(legend.position = "top")
save_plot(fc_plot, "GSE103940_A_vs_B_logFC_distribution", 8, 5.5)

scatter_plot <- ggplot(comparison, aes(logFC_A, logFC_B)) +
  geom_point(alpha = 0.4, size = 0.7, color = "#555555") +
  geom_abline(slope = 1, intercept = 0, color = "#C33C3C") +
  labs(x = "Analysis A logFC", y = "Analysis B logFC",
       title = sprintf("A/B logFC agreement: Pearson r = %.4f",
                       cor(comparison$logFC_A, comparison$logFC_B))) +
  theme_classic(base_size = 11)
save_plot(scatter_plot, "GSE103940_A_vs_B_logFC_scatter", 6.5, 6)

run_gsea <- function(tt, collection, gmt_file, method) {
  ranked <- tt[order(-tt$t, tt$gene), c("gene", "t")]
  stats <- ranked$t
  names(stats) <- ranked$gene
  pathways <- gmtPathways(gmt_file)
  result <- fgseaMultilevel(pathways = pathways, stats = stats, minSize = 10,
                            maxSize = 500, eps = 0, nproc = 1)
  result <- as.data.frame(result)
  result$leadingEdge <- vapply(result$leadingEdge, paste, collapse = ";",
                               FUN.VALUE = character(1))
  result$collection <- collection
  result$method <- method
  result$direction <- ifelse(result$NES > 0, "up_in_high_altitude", "down_in_high_altitude")
  result[order(result$padj, -abs(result$NES)),
         c("collection", "method", "pathway", "NES", "pval", "padj", "direction",
           "size", "ES", "log2err", "leadingEdge")]
}

gsea_results <- list()
gsea_consistency <- list()
for (collection in names(gmt_files)) {
  message("Running GSE103940 GSEA: ", collection, " analysis A")
  ga <- run_gsea(tt_a, collection, gmt_files[[collection]], "A_standard_TMM")
  message("Running GSE103940 GSEA: ", collection, " analysis B")
  gb <- run_gsea(tt_b, collection, gmt_files[[collection]], "B_non_globin_factor_TMM")
  write.csv(ga, file.path(out, paste0("GSE103940_analysis_A_", collection, "_GSEA.csv")),
            row.names = FALSE)
  write.csv(gb, file.path(out, paste0("GSE103940_analysis_B_", collection, "_GSEA.csv")),
            row.names = FALSE)
  joined <- merge(ga[, c("pathway", "NES", "pval", "padj")],
                  gb[, c("pathway", "NES", "pval", "padj")], by = "pathway",
                  suffixes = c("_A", "_B"))
  sig_a <- ga$pathway[ga$padj < 0.05]
  sig_b <- gb$pathway[gb$padj < 0.05]
  gsea_consistency[[collection]] <- data.frame(
    collection = collection, tested_pathways_A = nrow(ga), tested_pathways_B = nrow(gb),
    common_pathways = nrow(joined), NES_Pearson = cor(joined$NES_A, joined$NES_B),
    NES_Spearman = cor(joined$NES_A, joined$NES_B, method = "spearman"),
    significant_A = length(sig_a), significant_B = length(sig_b),
    significant_overlap = length(intersect(sig_a, sig_b)),
    significant_Jaccard = jaccard(sig_a, sig_b))
  write.csv(joined, file.path(out, paste0("GSE103940_A_vs_B_", collection,
                                          "_GSEA_comparison.csv")), row.names = FALSE)
  gsea_results[[paste0(collection, "_A")]] <- ga
  gsea_results[[paste0(collection, "_B")]] <- gb
}
gsea_consistency <- do.call(rbind, gsea_consistency)
write.csv(gsea_consistency, file.path(out, "GSE103940_A_vs_B_GSEA_consistency.csv"),
          row.names = FALSE)

metric_value <- setNames(metrics$value, metrics$metric)
freeze_pass <- metric_value[["logFC_Pearson"]] >= 0.98 &&
  metric_value[["moderated_t_Pearson"]] >= 0.98 &&
  metric_value[["primary_Jaccard"]] >= 0.80 &&
  metric_value[["direction_concordance"]] >= 0.98
freeze_decision <- data.frame(
  criterion = c("logFC Pearson >= 0.98", "moderated t Pearson >= 0.98",
                "primary DEG Jaccard >= 0.80", "direction concordance >= 0.98",
                "all freeze criteria"),
  observed = c(metric_value[["logFC_Pearson"]], metric_value[["moderated_t_Pearson"]],
               metric_value[["primary_Jaccard"]], metric_value[["direction_concordance"]],
               as.numeric(freeze_pass)),
  pass = c(metric_value[["logFC_Pearson"]] >= 0.98,
           metric_value[["moderated_t_Pearson"]] >= 0.98,
           metric_value[["primary_Jaccard"]] >= 0.80,
           metric_value[["direction_concordance"]] >= 0.98, freeze_pass)
)
write.csv(freeze_decision, file.path(out, "GSE103940_freeze_decision.csv"), row.names = FALSE)
writeLines(if (freeze_pass) {
  "PASS: Analysis A satisfies all prespecified globin-sensitivity criteria and may be frozen as the formal GSE103940 result."
} else {
  "FAIL: Analysis A does not satisfy all prespecified globin-sensitivity criteria and must not be frozen without further methodological review."
}, file.path(out, "GSE103940_freeze_decision.txt"))

config <- data.frame(
  item = c("input", "tested_gene_filter", "analysis_A", "analysis_B", "design", "contrast",
           "primary_threshold", "globin_genes_excluded_from_B_factor_estimation",
           "GSEA", "seed"),
  value = c(count_file, "edgeR filterByExpr using paired design; same genes in A and B",
            "standard TMM -> voom -> paired limma -> robust eBayes",
            "TMM factors estimated without major globins; applied to full matrix -> same voom/limma",
            "~ subject_id + condition", "high_altitude - plain",
            "BH-FDR < 0.05 and |log2FC| > 1.5", paste(globin_in_matrix, collapse = ";"),
            "MSigDB 2025.1.Hs Hallmark/Reactome/GO BP; fgseaMultilevel eps=0", "4103940")
)
write.csv(config, file.path(out, "GSE103940_final_analysis_config.csv"), row.names = FALSE)
input_paths <- c(count_file, meta_file, map_file, target_file, phase3_file, gmt_files)
write.csv(data.frame(path = input_paths, md5 = unname(tools::md5sum(input_paths))),
          file.path(out, "GSE103940_final_input_md5.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "GSE103940_final_sessionInfo.txt"))
saveRDS(list(design = design, metadata = meta, y_A = y_a, y_B = y_b,
             voom_A = v_a, voom_B = v_b, fit_A = fit_a, fit_B = fit_b,
             globin_genes = globin_in_matrix),
        file.path(out, "GSE103940_final_A_B_fit.rds"))

message("GSE103940 A/B complete. A primary=", length(primary_a),
        "; B primary=", length(primary_b), "; overlap=", length(intersect(primary_a, primary_b)),
        "; logFC r=", round(metric_value[["logFC_Pearson"]], 5),
        "; freeze=", freeze_pass)
