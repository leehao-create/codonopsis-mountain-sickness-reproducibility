#!/usr/bin/env Rscript

# Reviewer-requested MMP9 diagnostics from frozen 12-module assignments.
# This script does not rebuild or redefine the WGCNA modules.

args <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
input_file <- file.path(repo_root, "analysis", "WGCNA", "data", "processed",
                        "GSE75665_WGCNA_input_FPKM.csv")
assignment_file <- file.path(repo_root, "analysis", "WGCNA", "manuscript_primary", "results",
                             "WGCNA_module_assignments.csv")
eigengene_file <- file.path(repo_root, "analysis", "WGCNA", "manuscript_primary", "results",
                            "WGCNA_module_eigengenes.csv")
trait_file <- file.path(repo_root, "analysis", "WGCNA", "reviewer_checks", "results",
                        "WGCNA_binary_AMS_trait.csv")
output_file <- file.path(repo_root, "analysis", "WGCNA", "reviewer_checks", "reproduced_results",
                         "WGCNA_MMP9_diagnostics.csv")
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(WGCNA))

fpkm <- read.csv(input_file, row.names = 1, check.names = FALSE)
fpkm <- na.omit(fpkm)
fpkm <- fpkm[rowSums(fpkm) != 0, , drop = FALSE]
fpkm <- fpkm[apply(fpkm, 1, mad) > 1.5, , drop = FALSE]
expr <- as.data.frame(t(fpkm))

assignment <- read.csv(assignment_file, check.names = FALSE)
assignment <- assignment[match(colnames(expr), assignment$gene), ]
stopifnot(!anyNA(assignment$gene), identical(assignment$gene, colnames(expr)))

MEs <- read.csv(eigengene_file, row.names = 1, check.names = FALSE)
MEs <- MEs[rownames(expr), , drop = FALSE]
trait <- read.csv(trait_file, row.names = 1, check.names = FALSE)
trait <- trait[rownames(expr), , drop = FALSE]
stopifnot(identical(rownames(MEs), rownames(expr)), identical(rownames(trait), rownames(expr)))

gene <- "MMP9"
module <- assignment$module_color[assignment$gene == gene]
module_genes <- assignment$gene[assignment$module_color == module]
module_ME <- paste0("ME", module)

kME <- cor(expr[, gene], MEs[, module_ME], use = "pairwise.complete.obs", method = "pearson")
GS <- cor(expr[, gene], trait$AMS_status, use = "pairwise.complete.obs", method = "pearson")
GS_P <- corPvalueStudent(GS, nrow(expr))

signed_adjacency <- adjacency(expr[, module_genes, drop = FALSE],
                              power = 16, type = "signed", corFnc = "cor")
kWithin <- rowSums(signed_adjacency) - diag(signed_adjacency)
rank_within <- rank(-kWithin, ties.method = "min")

result <- data.frame(
  gene = gene,
  module = module,
  kME = unname(kME),
  traditional_GS = unname(GS),
  GS_P = unname(GS_P),
  kWithin_signed_adjacency = unname(kWithin[gene]),
  within_module_rank = unname(rank_within[gene]),
  module_size = length(module_genes),
  hub_interpretation = "Exploratory module context only; not a high-connectivity hub"
)
write.csv(result, output_file, row.names = FALSE)
cat("Wrote", output_file, "\n")
