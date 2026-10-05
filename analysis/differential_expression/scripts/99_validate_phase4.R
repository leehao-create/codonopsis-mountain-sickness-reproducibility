options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
module_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

g103 <- file.path(module_dir, "results", "GSE103940")
g756 <- file.path(module_dir, "results", "GSE75665")
g260 <- file.path(module_dir, "results", "GSE260910_confounding")

required <- c(
  file.path(g103, "GSE103940_analysis_A_DEG_all.csv"),
  file.path(g103, "GSE103940_analysis_A_DEG_significant.csv"),
  file.path(g103, "GSE103940_final_design_matrix.csv"),
  file.path(g103, "GSE103940_final_analysis_config.csv"),
  file.path(g756, "GSE75665_AMS_interaction_fixed_ranked_statistics.csv"),
  file.path(g756, "GSE75665_AMS_interaction_GSEA_all_collections.csv"),
  file.path(g756, "GSE75665_AMS_interaction_GSEA_FDR_significant_all_collections.csv"),
  file.path(g260, "GSE260910_sample_disease_batch_table.csv"),
  file.path(g260, "GSE260910_batch_confounding_report.md")
)
stopifnot(all(file.exists(required)), all(file.info(required)$size > 0))

a <- read.csv(file.path(g103, "GSE103940_analysis_A_DEG_all.csv"), check.names = FALSE)
sig <- read.csv(file.path(g103, "GSE103940_analysis_A_DEG_significant.csv"), check.names = FALSE)
interaction <- read.csv(file.path(module_dir, "data", "dependencies",
                                  "GSE75665_NCBI_counts_AMS_interaction_tested_genes.csv"),
                        check.names = FALSE)
stopifnot(nrow(a) == 14912L, nrow(sig) == 296L,
          sum(sig$logFC > 0) == 283L, sum(sig$logFC < 0) == 13L,
          sum(interaction$adj.P.Val < 0.05, na.rm = TRUE) == 0L)

cat("PASS: frozen bulk release tables are present and internally consistent\n")
