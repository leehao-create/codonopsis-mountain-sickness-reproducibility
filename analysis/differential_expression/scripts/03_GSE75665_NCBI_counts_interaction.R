options(stringsAsFactors = FALSE)
set.seed(75665)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
module_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(ggplot2)
  library(pheatmap)
  library(fgsea)
})

out <- Sys.getenv("DE_OUTPUT_DIR", file.path(module_dir, "reproduced_results", "GSE75665_interaction"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
count_file <- file.path(module_dir, "data", "processed", "GSE75665_raw_counts_GRCh38.p13_NCBI.tsv.gz")
meta_file <- file.path(module_dir, "data", "metadata", "GSE75665_subject_pairing.csv")
map_file <- Sys.getenv(
  "NCBI_GENE_MAP",
  file.path(module_dir, "data", "reference", "NCBI_GeneID_to_symbol.txt")
)
target_file <- file.path(module_dir, "..", "network_pharmacology", "data", "derived",
                         "Codonopsis_mountain_sickness_54_targets.csv")
gmt_file <- file.path(Sys.getenv("MSIGDB_DIR", file.path(module_dir, "data", "gene_sets")),
                      "h.all.v2025.1.Hs.symbols.gmt")
rpkm_results_file <- file.path(module_dir, "data", "dependencies",
                               "GSE75665_AMS_interaction_tested_genes.csv")
if (!file.exists(map_file)) {
  stop("NCBI gene mapping is not redistributed. Run scripts/download_ncbi_gene_mapping.sh ",
       "or set NCBI_GENE_MAP; see the module README.")
}
stopifnot(dir.exists(out), file.exists(count_file), file.exists(meta_file),
          file.exists(target_file), file.exists(gmt_file),
          file.exists(rpkm_results_file))

save_plot <- function(plot, stem, width = 8, height = 6) {
  ggsave(file.path(out, paste0(stem, ".pdf")), plot, width = width, height = height)
  ggsave(file.path(out, paste0(stem, ".png")), plot, width = width, height = height,
         dpi = 300)
}

meta <- read.csv(meta_file, check.names = FALSE)
stopifnot(nrow(meta) == 20L, !anyDuplicated(meta$gsm), all(table(meta$subject_id) == 2L))
count_df <- read.delim(count_file, check.names = FALSE)
stopifnot(colnames(count_df)[1] == "GeneID", identical(colnames(count_df)[-1], meta$gsm),
          !anyDuplicated(count_df$GeneID))
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

subject_factor <- factor(meta$subject_id)
altitude_high <- as.integer(meta$altitude_status == "high_altitude")
ams_indicator <- as.integer(meta$ams_status == "AMS")
design <- cbind(model.matrix(~ 0 + subject_factor), altitude_high = altitude_high,
                AMS_x_altitude = ams_indicator * altitude_high)
colnames(design)[seq_len(nlevels(subject_factor))] <- paste0("subject_", levels(subject_factor))
rownames(design) <- meta$gsm
stopifnot(qr(design)$rank == ncol(design), ncol(design) == 12L)

y_all <- DGEList(counts = counts_symbol)
four_groups <- interaction(meta$ams_status, meta$altitude_status, drop = TRUE)
keep <- filterByExpr(y_all, group = four_groups)
y <- y_all[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
v <- voom(y, design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit, trend = FALSE, robust = TRUE)
tt <- topTable(fit, coef = "AMS_x_altitude", number = Inf, adjust.method = "BH", sort.by = "P")
tt$gene <- rownames(tt)

subjects <- unique(meta$subject_id)
delta <- sapply(subjects, function(subject) {
  idx <- which(meta$subject_id == subject)
  h <- idx[meta$altitude_status[idx] == "high_altitude"]
  p <- idx[meta$altitude_status[idx] == "plain"]
  v$E[, h] - v$E[, p]
})
colnames(delta) <- subjects
subject_status <- vapply(subjects, function(subject) unique(meta$ams_status[meta$subject_id == subject]),
                         character(1))
direct_interaction <- rowMeans(delta[, subject_status == "AMS", drop = FALSE]) -
  rowMeans(delta[, subject_status == "non_AMS", drop = FALSE])
unweighted_fit <- lmFit(v$E, design)
unweighted_interaction <- unweighted_fit$coefficients[, "AMS_x_altitude"]
stopifnot(max(abs(unweighted_interaction - direct_interaction[names(unweighted_interaction)])) < 1e-10)

tt$mean_delta_AMS <- rowMeans(delta[tt$gene, subject_status == "AMS", drop = FALSE])
tt$mean_delta_non_AMS <- rowMeans(delta[tt$gene, subject_status == "non_AMS", drop = FALSE])
tt$direction <- ifelse(tt$logFC > 0, "larger_altitude_response_in_AMS",
                       ifelse(tt$logFC < 0, "smaller_altitude_response_in_AMS", "no_difference"))
tt$passes_FDR_0_05 <- tt$adj.P.Val < 0.05
tt$passes_primary_threshold <- tt$adj.P.Val < 0.05 & abs(tt$logFC) > 1.5
tt <- tt[, c("gene", setdiff(colnames(tt), "gene"))]

all_genes <- data.frame(gene = rownames(counts_symbol), expression_eligible = keep,
                        tested = keep, stringsAsFactors = FALSE)
all_genes <- merge(all_genes, tt, by = "gene", all.x = TRUE, sort = FALSE)
all_genes <- all_genes[match(rownames(counts_symbol), all_genes$gene), ]
stopifnot(identical(all_genes$gene, rownames(counts_symbol)))
write.csv(all_genes, file.path(out, "GSE75665_NCBI_counts_AMS_interaction_all_genes.csv"),
          row.names = FALSE)
write.csv(tt, file.path(out, "GSE75665_NCBI_counts_AMS_interaction_tested_genes.csv"),
          row.names = FALSE)
write.csv(tt[tt$passes_primary_threshold, ],
          file.path(out, "GSE75665_NCBI_counts_AMS_interaction_DEG_significant.csv"),
          row.names = FALSE)
write.csv(tt[tt$passes_FDR_0_05, ],
          file.path(out, "GSE75665_NCBI_counts_AMS_interaction_FDR_only.csv"), row.names = FALSE)

counts <- data.frame(
  criterion = c("BH-FDR < 0.05 and |interaction log2FC| > 1.5", "BH-FDR < 0.05 only"),
  total = c(sum(tt$passes_primary_threshold), sum(tt$passes_FDR_0_05)),
  positive_interaction = c(sum(tt$passes_primary_threshold & tt$logFC > 0),
                           sum(tt$passes_FDR_0_05 & tt$logFC > 0)),
  negative_interaction = c(sum(tt$passes_primary_threshold & tt$logFC < 0),
                           sum(tt$passes_FDR_0_05 & tt$logFC < 0))
)
write.csv(counts, file.path(out, "GSE75665_NCBI_counts_AMS_interaction_DEG_counts.csv"),
          row.names = FALSE)
write.csv(data.frame(
  contrast = "(AMS_high - AMS_plain) - (nonAMS_high - nonAMS_plain)",
  design_rows = nrow(design), design_columns = ncol(design), design_rank = qr(design)$rank,
  maximum_absolute_unweighted_model_vs_direct_difference =
    max(abs(unweighted_interaction - direct_interaction[names(unweighted_interaction)])),
  voom_weighted_vs_direct_Pearson_correlation = cor(tt$logFC, direct_interaction[tt$gene]),
  maximum_absolute_voom_weighted_vs_direct_difference =
    max(abs(tt$logFC - direct_interaction[tt$gene]))
), file.path(out, "GSE75665_NCBI_counts_AMS_interaction_validation.csv"), row.names = FALSE)

sample_summary <- data.frame(
  sample_id = meta$gsm, subject_id = meta$subject_id, ams_status = meta$ams_status,
  altitude_status = meta$altitude_status, raw_library_size = y$samples$lib.size,
  TMM_norm_factor = y$samples$norm.factors,
  effective_library_size = y$samples$lib.size * y$samples$norm.factors,
  voom_logCPM_mean = colMeans(v$E), voom_logCPM_median = apply(v$E, 2, median)
)
write.csv(sample_summary, file.path(out, "GSE75665_NCBI_counts_sample_summary.csv"),
          row.names = FALSE)
write.csv(data.frame(rule = c("mapped unique gene symbols", "filterByExpr four groups"),
                     genes = c(nrow(counts_symbol), sum(keep))),
          file.path(out, "GSE75665_NCBI_counts_filter_summary.csv"), row.names = FALSE)

targets <- unique(trimws(read.csv(target_file, check.names = FALSE)[[1]]))
stopifnot(length(targets) == 54L)
target_results <- data.frame(gene = targets, in_input = targets %in% rownames(counts_symbol),
                             expression_eligible = targets %in% rownames(v$E))
target_results <- merge(target_results, tt, by = "gene", all.x = TRUE, sort = FALSE)
target_results <- target_results[match(targets, target_results$gene), ]
write.csv(target_results,
          file.path(out, "GSE75665_NCBI_counts_Codonopsis_targets_interaction_results.csv"),
          row.names = FALSE)
focus <- c("MMP9", "NR3C2", "HBB", "SNCA", "ABCB1", "BCL2", "CA2")
write.csv(target_results[match(focus, target_results$gene), ],
          file.path(out, "GSE75665_NCBI_counts_focus_candidates_interaction_results.csv"),
          row.names = FALSE)

plot_df <- tt
plot_df$category <- "not_primary"
plot_df$category[plot_df$passes_primary_threshold & plot_df$logFC > 0] <- "positive_interaction"
plot_df$category[plot_df$passes_primary_threshold & plot_df$logFC < 0] <- "negative_interaction"
plot_df$minus_log10_FDR <- -log10(pmax(plot_df$adj.P.Val, .Machine$double.xmin))
volcano <- ggplot(plot_df, aes(logFC, minus_log10_FDR, color = category)) +
  geom_point(alpha = 0.65, size = 1.15) +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = 2, color = "grey35") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey35") +
  scale_color_manual(values = c(negative_interaction = "#2878B5", not_primary = "#A7A9AC",
                                positive_interaction = "#C33C3C")) +
  labs(x = "Interaction log2FC", y = "-log10(BH-FDR)", color = NULL,
       title = "GSE75665 count-based AMS-specific altitude response") +
  theme_classic(base_size = 11) + theme(legend.position = "top")
save_plot(volcano, "GSE75665_NCBI_counts_AMS_interaction_volcano", 8, 6)

heat_genes <- tt$gene[tt$passes_primary_threshold]
heatmap_basis <- "primary significant interaction genes"
if (length(heat_genes) < 2L) {
  heat_genes <- head(tt$gene, 50L)
  heatmap_basis <- "top 50 genes by raw interaction P value; not a significance claim"
}
heat <- delta[heat_genes, , drop = FALSE]
heat_z <- t(scale(t(heat)))
heat_z[!is.finite(heat_z)] <- 0
heat_ann <- data.frame(AMS_status = factor(subject_status, levels = c("non_AMS", "AMS")),
                       row.names = subjects)
write.csv(data.frame(gene = heat_genes, heatmap_basis = heatmap_basis),
          file.path(out, "GSE75665_NCBI_counts_AMS_interaction_heatmap_genes.csv"),
          row.names = FALSE)
pdf(file.path(out, "GSE75665_NCBI_counts_AMS_interaction_top_genes_heatmap.pdf"),
    width = 8, height = 10)
pheatmap(heat_z, annotation_col = heat_ann, cluster_cols = FALSE, fontsize_row = 6,
         main = paste0("High - plain voom change (row z-score)\n", heatmap_basis))
dev.off()
png(file.path(out, "GSE75665_NCBI_counts_AMS_interaction_top_genes_heatmap.png"),
    width = 2400, height = 3000, res = 300)
pheatmap(heat_z, annotation_col = heat_ann, cluster_cols = FALSE, fontsize_row = 6,
         main = paste0("High - plain voom change (row z-score)\n", heatmap_basis))
dev.off()

ranked <- tt[, c("gene", "t", "logFC", "P.Value", "adj.P.Val")]
ranked <- ranked[order(ranked$t, decreasing = TRUE), ]
write.csv(ranked, file.path(out, "GSE75665_NCBI_counts_AMS_interaction_ranked_genes.csv"),
          row.names = FALSE)
pathways <- gmtPathways(gmt_file)
stats <- ranked$t
names(stats) <- ranked$gene
stats <- sort(stats, decreasing = TRUE)
gsea <- fgseaMultilevel(pathways = pathways, stats = stats, minSize = 10, maxSize = 500,
                        eps = 0)
gsea <- as.data.frame(gsea)
gsea$leadingEdge <- vapply(gsea$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
gsea <- gsea[order(gsea$padj, -abs(gsea$NES)), ]
write.csv(gsea, file.path(out, "GSE75665_NCBI_counts_AMS_interaction_Hallmark_GSEA.csv"),
          row.names = FALSE)
write.csv(gsea[gsea$padj < 0.05, ],
          file.path(out, "GSE75665_NCBI_counts_AMS_interaction_Hallmark_GSEA_significant.csv"),
          row.names = FALSE)
top_gsea <- head(gsea, 15L)
top_gsea$pathway_label <- factor(gsub("HALLMARK_", "", top_gsea$pathway),
                                 levels = rev(gsub("HALLMARK_", "", top_gsea$pathway)))
gsea_plot <- ggplot(top_gsea, aes(NES, pathway_label, color = padj < 0.05)) +
  geom_point(aes(size = -log10(pmax(padj, .Machine$double.xmin)))) +
  geom_vline(xintercept = 0, color = "grey50") +
  scale_color_manual(values = c(`FALSE` = "#777777", `TRUE` = "#C33C3C")) +
  labs(x = "Normalized enrichment score", y = NULL, color = "FDR < 0.05",
       size = "-log10(FDR)", title = "Count-based Hallmark GSEA: AMS interaction") +
  theme_classic(base_size = 10)
save_plot(gsea_plot, "GSE75665_NCBI_counts_AMS_interaction_Hallmark_GSEA", 9, 6.5)

rpkm_tt <- read.csv(rpkm_results_file, check.names = FALSE)
common <- merge(tt[, c("gene", "logFC", "P.Value", "adj.P.Val", "passes_primary_threshold")],
                rpkm_tt[, c("gene", "logFC", "P.Value", "adj.P.Val", "passes_primary_threshold")],
                by = "gene", suffixes = c("_NCBI_counts", "_submitter_RPKM"))
common$same_direction <- sign(common$logFC_NCBI_counts) == sign(common$logFC_submitter_RPKM)
write.csv(common, file.path(out, "GSE75665_counts_vs_RPKM_interaction_comparison.csv"),
          row.names = FALSE)
comparison_summary <- data.frame(
  metric = c("common_tested_genes", "Pearson_logFC", "Spearman_logFC",
             "overall_sign_concordance", "count_primary", "RPKM_primary", "primary_overlap"),
  value = c(nrow(common), cor(common$logFC_NCBI_counts, common$logFC_submitter_RPKM),
            cor(common$logFC_NCBI_counts, common$logFC_submitter_RPKM, method = "spearman"),
            mean(common$same_direction), sum(common$passes_primary_threshold_NCBI_counts),
            sum(common$passes_primary_threshold_submitter_RPKM),
            sum(common$passes_primary_threshold_NCBI_counts &
                  common$passes_primary_threshold_submitter_RPKM))
)
write.csv(comparison_summary,
          file.path(out, "GSE75665_counts_vs_RPKM_interaction_summary.csv"), row.names = FALSE)

config <- data.frame(
  item = c("primary_input", "mapping", "filter", "normalization", "design", "contrast",
           "eBayes", "primary_threshold", "GSEA", "seed"),
  value = c(count_file, "NCBI Gene ID to symbol; duplicate symbols summed",
            "edgeR filterByExpr with four AMS-by-altitude groups", "TMM then voom",
            "10 subject fixed effects + altitude_high + AMS_x_altitude",
            "(AMS_high-AMS_plain)-(nonAMS_high-nonAMS_plain)", "robust=TRUE",
            "BH-FDR < 0.05 and |interaction log2FC| > 1.5",
            "MSigDB 2025.1.Hs Hallmark; moderated t; fgseaMultilevel eps=0", "75665")
)
write.csv(config, file.path(out, "GSE75665_NCBI_counts_analysis_config.csv"), row.names = FALSE)
write.csv(data.frame(path = c(count_file, meta_file, map_file, target_file, gmt_file),
                     md5 = unname(tools::md5sum(c(count_file, meta_file, map_file,
                                                  target_file, gmt_file)))),
          file.path(out, "GSE75665_NCBI_counts_input_md5.csv"), row.names = FALSE)
capture.output(sessionInfo(),
               file = file.path(out, "GSE75665_NCBI_counts_phase3_sessionInfo.txt"))
saveRDS(list(design = design, fit = fit, metadata = meta, voom = v, subject_delta = delta),
        file.path(out, "GSE75665_NCBI_counts_AMS_interaction_fit.rds"))

message("GSE75665 count interaction complete: ", counts$total[1], " primary; ",
        counts$total[2], " FDR-only; Hallmark FDR pathways=", sum(gsea$padj < 0.05))
