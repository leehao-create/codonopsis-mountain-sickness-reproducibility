#!/usr/bin/env Rscript

# Repository-portable transcription of the archived manuscript workflow.
# Core analytical calls are preserved, including the absence of an explicit
# networkType argument in blockwiseModules().

args <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args) >= 1L) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
out_dir <- if (length(args) >= 2L) args[[2]] else file.path(
  repo_root, "analysis", "WGCNA", "manuscript_primary", "reproduced_results"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(repo_root, "analysis", "WGCNA", "data", "processed",
                        "GSE75665_WGCNA_input_FPKM.csv")
metadata_file <- file.path(repo_root, "analysis", "WGCNA", "manuscript_primary", "data",
                           "GSE75665_original_4group_metadata.csv")

suppressPackageStartupMessages({
  library(WGCNA)
  library(stringr)
})
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

dataExpr <- read.csv(input_file, row.names = 1, check.names = FALSE)
dataExpr <- na.omit(dataExpr)
dataExpr <- dataExpr[rowSums(dataExpr) != 0, , drop = FALSE]
stopifnot(nrow(dataExpr) == 20751L, ncol(dataExpr) == 20L)

m.mad <- apply(dataExpr, 1, mad)
dataExprVar <- dataExpr[m.mad > 1.5, , drop = FALSE]
stopifnot(nrow(dataExprVar) == 2810L)
dataExpr <- as.data.frame(t(dataExprVar))

gsg <- goodSamplesGenes(dataExpr, verbose = 3)
if (!gsg$allOK) {
  dataExpr <- dataExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
}
stopifnot(nrow(dataExpr) == 20L, ncol(dataExpr) == 2810L)

sampleTree <- hclust(dist(dataExpr), method = "average")
pdf(file.path(out_dir, "WGCNA_sample_clustering.pdf"), width = 10, height = 6)
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")
dev.off()

powers <- c(1:10, seq(from = 12, to = 30, by = 2))
sft <- pickSoftThreshold(dataExpr, powerVector = powers,
                         networkType = "signed", verbose = 5)
write.csv(sft$fitIndices,
          file.path(out_dir, "WGCNA_soft_threshold_fitIndices.csv"),
          row.names = FALSE)

nGenes <- ncol(dataExpr)
nSamples <- nrow(dataExpr)
net <- blockwiseModules(
  dataExpr,
  power = 16,
  maxBlockSize = nGenes,
  TOMType = "signed",
  minModuleSize = 20,
  reassignThreshold = 0,
  mergeCutHeight = 0.05,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  corType = "pearson",
  loadTOMs = TRUE,
  verbose = 3
)

moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)
stopifnot(length(unique(moduleLabels)) == 12L)

MEs <- net$MEs
MEs_col <- MEs
colnames(MEs_col) <- paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs), "ME", ""))
))
MEs_col <- orderMEs(MEs_col)

metadata <- read.csv(metadata_file, check.names = FALSE)
sample_col <- if ("sample_id" %in% names(metadata)) "sample_id" else names(metadata)[[1]]
group_col <- if ("group" %in% names(metadata)) "group" else stop("Metadata lacks the original four-level group column")
rownames(metadata) <- metadata[[sample_col]]
metadata <- metadata[rownames(dataExpr), , drop = FALSE]
stopifnot(identical(rownames(metadata), rownames(dataExpr)))

design <- model.matrix(~ 0 + factor(metadata[[group_col]]))
colnames(design) <- levels(factor(metadata[[group_col]]))
rownames(design) <- rownames(dataExpr)
modTraitCor <- cor(MEs_col, design, use = "p")
modTraitP <- corPvalueStudent(modTraitCor, nSamples)

assignments <- data.frame(
  gene = colnames(dataExpr), module_label = moduleLabels,
  module_color = moduleColors, stringsAsFactors = FALSE
)
write.csv(assignments, file.path(out_dir, "WGCNA_module_assignments.csv"), row.names = FALSE)
write.csv(MEs_col, file.path(out_dir, "WGCNA_module_eigengenes.csv"))
write.csv(modTraitCor, file.path(out_dir, "WGCNA_module_trait_4group_correlations.csv"))
write.csv(modTraitP, file.path(out_dir, "WGCNA_module_trait_4group_pvalues.csv"))
saveRDS(net, file.path(out_dir, "WGCNA_network_object.rds"))
capture.output(sessionInfo(), file = file.path(out_dir, "WGCNA_sessionInfo.txt"))

cat("Completed manuscript-primary workflow:", nGenes, "genes,",
    length(unique(moduleLabels)), "modules. Outputs:", out_dir, "\n")
