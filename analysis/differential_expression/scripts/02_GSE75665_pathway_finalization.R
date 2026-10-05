options(stringsAsFactors = FALSE)
set.seed(75665)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
module_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
suppressPackageStartupMessages({
  library(fgsea)
  library(ggplot2)
})

out <- Sys.getenv("DE_OUTPUT_DIR", file.path(module_dir, "reproduced_results", "GSE75665_GSEA"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rank_file <- file.path(module_dir, "data", "dependencies",
                       "GSE75665_NCBI_counts_AMS_interaction_tested_genes.csv")
msigdb_dir <- Sys.getenv("MSIGDB_DIR", file.path(module_dir, "data", "gene_sets"))
gmt_files <- c(
  Hallmark = file.path(msigdb_dir, "h.all.v2025.1.Hs.symbols.gmt"),
  Reactome = file.path(msigdb_dir, "c2.cp.reactome.v2025.1.Hs.symbols.gmt"),
  GO_BP = file.path(msigdb_dir, "c5.go.bp.v2025.1.Hs.symbols.gmt")
)
stopifnot(dir.exists(out), file.exists(rank_file), all(file.exists(gmt_files)))

rank_table <- read.csv(rank_file, check.names = FALSE)
required <- c("gene", "logFC", "t", "P.Value", "adj.P.Val",
              "passes_FDR_0_05", "passes_primary_threshold")
stopifnot(all(required %in% colnames(rank_table)), !anyDuplicated(rank_table$gene),
          nrow(rank_table) == 13742L, all(is.finite(rank_table$t)),
          sum(rank_table$adj.P.Val < 0.05, na.rm = TRUE) == 0L,
          sum(rank_table$passes_FDR_0_05, na.rm = TRUE) == 0L,
          sum(rank_table$passes_primary_threshold, na.rm = TRUE) == 0L)

ranked <- rank_table[order(-rank_table$t, rank_table$gene), c("gene", "t")]
stats <- ranked$t
names(stats) <- ranked$gene
stopifnot(!anyDuplicated(names(stats)))

run_gsea <- function(collection, gmt_file) {
  pathways <- gmtPathways(gmt_file)
  result <- fgseaMultilevel(pathways = pathways, stats = stats, minSize = 10,
                            maxSize = 500, eps = 0, nproc = 1)
  result <- as.data.frame(result)
  result$leadingEdge <- vapply(result$leadingEdge, paste, collapse = ";",
                               FUN.VALUE = character(1))
  result$collection <- collection
  result$direction <- ifelse(
    result$NES > 0,
    "larger_altitude_response_in_AMS",
    "smaller_or_reversed_altitude_response_in_AMS"
  )
  result <- result[order(result$padj, -abs(result$NES), result$pathway), ]
  result[, c("collection", "pathway", "NES", "pval", "padj", "direction",
             "size", "ES", "log2err", "leadingEdge")]
}

all_results <- list()
summary_rows <- list()
for (collection in names(gmt_files)) {
  message("Running fixed-rank GSE75665 GSEA: ", collection)
  result <- run_gsea(collection, gmt_files[[collection]])
  significant <- result[!is.na(result$padj) & result$padj < 0.05, ]
  write.csv(result,
            file.path(out, paste0("GSE75665_AMS_interaction_", collection, "_GSEA_all.csv")),
            row.names = FALSE)
  write.csv(significant,
            file.path(out, paste0("GSE75665_AMS_interaction_", collection,
                                  "_GSEA_FDR_significant.csv")),
            row.names = FALSE)
  all_results[[collection]] <- result
  summary_rows[[collection]] <- data.frame(
    collection = collection,
    tested_pathways = nrow(result),
    FDR_significant = nrow(significant),
    positive_NES_FDR_significant = sum(significant$NES > 0),
    negative_NES_FDR_significant = sum(significant$NES < 0),
    minimum_FDR = if (nrow(result)) min(result$padj, na.rm = TRUE) else NA_real_
  )

  plot_data <- head(significant, 20L)
  if (!nrow(plot_data)) {
    plot_data <- head(result, 20L)
    plot_subtitle <- "No pathways reached BH-FDR < 0.05; top ranked pathways shown"
  } else {
    plot_subtitle <- "BH-FDR < 0.05 pathways (top 20 if applicable)"
  }
  plot_data$label <- factor(gsub("_", " ", plot_data$pathway),
                            levels = rev(gsub("_", " ", plot_data$pathway)))
  p <- ggplot(plot_data, aes(NES, label, size = -log10(pmax(padj, .Machine$double.xmin)),
                             color = direction)) +
    geom_point(alpha = 0.85) +
    geom_vline(xintercept = 0, color = "#777777", linewidth = 0.35) +
    scale_color_manual(
      values = c(larger_altitude_response_in_AMS = "#B33A3A",
                 smaller_or_reversed_altitude_response_in_AMS = "#256E9B"),
      labels = c(larger_altitude_response_in_AMS = "Larger response in AMS",
                 smaller_or_reversed_altitude_response_in_AMS = "Smaller/reversed in AMS")
    ) +
    labs(x = "Normalized enrichment score (NES)", y = NULL,
         size = expression(-log[10](FDR)), color = NULL,
         title = paste("GSE75665 AMS interaction GSEA:", collection),
         subtitle = plot_subtitle) +
    theme_classic(base_size = 10) +
    theme(legend.position = "bottom")
  ggsave(file.path(out, paste0("GSE75665_AMS_interaction_", collection, "_GSEA.pdf")),
         p, width = 10, height = 7)
  ggsave(file.path(out, paste0("GSE75665_AMS_interaction_", collection, "_GSEA.png")),
         p, width = 10, height = 7, dpi = 300)
}

combined <- do.call(rbind, all_results)
rownames(combined) <- NULL
combined_sig <- combined[!is.na(combined$padj) & combined$padj < 0.05, ]
summary_table <- do.call(rbind, summary_rows)
rownames(summary_table) <- NULL
write.csv(combined, file.path(out, "GSE75665_AMS_interaction_GSEA_all_collections.csv"),
          row.names = FALSE)
write.csv(combined_sig,
          file.path(out, "GSE75665_AMS_interaction_GSEA_FDR_significant_all_collections.csv"),
          row.names = FALSE)
write.csv(summary_table, file.path(out, "GSE75665_AMS_interaction_GSEA_summary.csv"),
          row.names = FALSE)

priority_hallmark <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_HEME_METABOLISM"
)
hallmark <- all_results[["Hallmark"]]
priority <- hallmark[match(priority_hallmark, hallmark$pathway), ]
stopifnot(identical(priority$pathway, priority_hallmark))
priority$FDR_significant <- !is.na(priority$padj) & priority$padj < 0.05
write.csv(priority, file.path(out, "GSE75665_priority_Hallmark_pathways.csv"),
          row.names = FALSE)

config <- data.frame(
  item = c("rank_input", "biological_contrast", "ranking_statistic", "gene_level_result",
           "collections", "MSigDB_release", "algorithm", "gene_set_size",
           "pathway_significance", "seed"),
  value = c(rank_file,
            "(AMS_high - AMS_plain) - (nonAMS_high - nonAMS_plain)",
            "moderated t statistic from the fixed Phase 3 subject-aware interaction model",
            "0 genes at BH-FDR < 0.05; pathway analysis does not alter this conclusion",
            "Hallmark;Reactome;GO Biological Process", "2025.1.Hs",
            "fgseaMultilevel; eps=0; nproc=1", "minSize=10; maxSize=500",
            "BH-FDR < 0.05 within each prespecified collection", "75665")
)
write.csv(config, file.path(out, "GSE75665_AMS_interaction_GSEA_config.csv"),
          row.names = FALSE)
input_paths <- c(rank_file, gmt_files)
write.csv(data.frame(path = input_paths, md5 = unname(tools::md5sum(input_paths))),
          file.path(out, "GSE75665_AMS_interaction_GSEA_input_md5.csv"), row.names = FALSE)
write.csv(ranked, file.path(out, "GSE75665_AMS_interaction_fixed_ranked_statistics.csv"),
          row.names = FALSE)
capture.output(sessionInfo(),
               file = file.path(out, "GSE75665_AMS_interaction_GSEA_sessionInfo.txt"))

message("GSE75665 pathway finalization complete. Significant pathways: ",
        nrow(combined_sig), " across three separately adjusted collections.")
