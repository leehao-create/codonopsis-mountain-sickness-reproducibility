options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)

suppressPackageStartupMessages(library(WGCNA))

out_dir <- file.path(repo_root, "analysis", "WGCNA", "reviewer_checks", "results")
cor_file <- file.path(out_dir, "WGCNA_module_AMS_correlation_matrix.csv")
p_file <- file.path(out_dir, "WGCNA_module_AMS_pvalue_matrix.csv")

module_cor <- as.matrix(read.csv(cor_file, row.names = 1, check.names = FALSE))
module_p <- as.matrix(read.csv(p_file, row.names = 1, check.names = FALSE))

stopifnot(identical(dim(module_cor), dim(module_p)))
stopifnot(identical(rownames(module_cor), rownames(module_p)))
stopifnot(ncol(module_cor) == 1L, colnames(module_cor) == "AMS_status")

module_fdr <- p.adjust(module_p[, 1], method = "BH")
stopifnot(all(module_fdr >= 0.05))

cell_text <- matrix(
  sprintf("r = %.2f\nP = %.3f", module_cor[, 1], module_p[, 1]),
  nrow = nrow(module_cor), ncol = 1,
  dimnames = dimnames(module_cor)
)

draw_heatmap <- function() {
  par(mar = c(6.4, 10.5, 4.2, 3.8), xpd = NA, family = "Arial")
  labeledHeatmap(
    Matrix = module_cor,
    xLabels = "AMS status",
    yLabels = rownames(module_cor),
    ySymbols = rownames(module_cor),
    colorLabels = FALSE,
    colors = blueWhiteRed(101),
    textMatrix = cell_text,
    setStdMargins = FALSE,
    cex.text = 0.76,
    cex.lab = 0.82,
    xLabelsAngle = 0,
    xLabelsAdj = 0.5,
    zlim = c(-1, 1),
    legendLabel = "Pearson r",
    cex.legendLabel = 0.80,
    main = "WGCNA module-AMS status correlations"
  )
  mtext("C", side = 3, line = 2.65, adj = -0.24, cex = 1.15, font = 2)
  mtext(
    "Binary trait: AMS = 1, non-AMS = 0",
    side = 1, line = 2.8, cex = 0.78, font = 2
  )
  mtext(
    "None of the modules remained significant after BH correction.",
    side = 1, line = 4.0, cex = 0.74
  )
}

figure_dir <- file.path(repo_root, "analysis", "WGCNA", "reviewer_checks", "reproduced_figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
pdf_file <- file.path(figure_dir, "Figure3C_WGCNA_AMS_binary_trait_heatmap.pdf")
png_file <- file.path(figure_dir, "Figure3C_WGCNA_AMS_binary_trait_heatmap.png")

cairo_pdf(pdf_file, width = 5.4, height = 7.1, family = "Arial")
draw_heatmap()
dev.off()

png(png_file, width = 5.4, height = 7.1, units = "in", res = 600,
    type = "cairo", family = "Arial", bg = "white")
draw_heatmap()
dev.off()

cat("PDF:", pdf_file, "\n")
cat("PNG:", png_file, "\n")
cat("Minimum BH-FDR:", format(min(module_fdr), digits = 6), "\n")
