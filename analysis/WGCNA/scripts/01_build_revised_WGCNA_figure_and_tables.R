options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
module_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
source_dir <- Sys.getenv("WGCNA_RESULT_DIR", file.path(module_dir, "results"))
output_dir <- Sys.getenv("WGCNA_FIGURE_DIR", file.path(module_dir, "reproduced_figures"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_result <- function(name, ...) {
  read.csv(file.path(source_dir, name), check.names = FALSE, ...)
}

soft_main <- read_result("WGCNA_soft_threshold_candidates.csv")
# The archived candidate table already contains powers 1-30 plus the 31-50
# diagnostic extension, with `within_prespecified_range` marking their roles.
soft <- soft_main[order(soft_main$Power), ]
stopifnot(identical(soft$Power, 1:50))
soft$signed_fit <- -sign(soft$slope) * soft$SFT.R.sq

module_sizes <- read_result("WGCNA_module_sizes.csv")
module_sizes <- module_sizes[module_sizes$module_colors != "grey", ]
MEs <- read_result("WGCNA_module_eigengenes.csv", row.names = 1)
mixed <- read_result("WGCNA_module_trait_mixed_effects.csv")
network_object <- readRDS(file.path(source_dir, "WGCNA_primary_network.rds"))

module_order <- c("turquoise", "blue", "brown")
module_palette <- c(turquoise = "#42B8B0", blue = "#3B6FB6", brown = "#9A6546", grey = "#BDBDBD")
effect_order <- c(
  "AMS_statusAMS",
  "altitude_statushigh_altitude",
  "AMS_statusAMS:altitude_statushigh_altitude"
)
effect_labels <- c("AMS at plain", "Altitude in non-AMS", "AMS x altitude")

diverging_color <- function(x, limit) {
  pal <- grDevices::colorRampPalette(c("#7899C2", "#FAFAFA", "#D98679"))(201)
  idx <- round((pmax(-limit, pmin(limit, x)) + limit) / (2 * limit) * 200) + 1
  pal[idx]
}

panel_label <- function(label) {
  mtext(label, side = 3, line = 1.15, adj = -0.075, font = 2, cex = 1.25)
}

draw_soft_threshold <- function() {
  old <- par(mar = c(4.4, 4.6, 3.0, 4.6), xaxs = "i")
  on.exit(par(old))
  main_idx <- soft$Power <= 30
  plot(soft$Power, soft$signed_fit, type = "n", xlim = c(1, 50), ylim = c(-0.65, 1),
       xlab = "Soft-threshold power", ylab = expression("Signed scale-free fit " * R^2),
       main = "Soft-threshold diagnostics", cex.main = 1.05)
  abline(h = 0.85, col = "#888888", lty = 2, lwd = 1)
  abline(v = 30, col = "#B22222", lty = 3, lwd = 1)
  lines(soft$Power[main_idx], soft$signed_fit[main_idx], col = "#2C7BB6", lwd = 1.4)
  points(soft$Power[main_idx], soft$signed_fit[main_idx], pch = 16, col = "#2C7BB6", cex = 0.52)
  lines(soft$Power[!main_idx], soft$signed_fit[!main_idx], col = "#8C8C8C", lwd = 1.1)
  points(soft$Power[!main_idx], soft$signed_fit[!main_idx], pch = 1, col = "#777777", cex = 0.48)
  points(30, soft$signed_fit[soft$Power == 30], pch = 21, bg = "#D95F4C", col = "#8B1A1A", cex = 1.0)
  text(31.0, 0.69, expression("Selected " * beta * " = 30"),
       adj = c(0, 0.5), cex = 0.76, col = "#8B1A1A")
  text(31.0, 0.58, expression(R^2 * " = 0.103"),
       adj = c(0, 0.5), cex = 0.76, col = "#8B1A1A")
  text(31.0, 0.49, "criterion not reached",
       adj = c(0, 0.5), cex = 0.72, col = "#8B1A1A")
  par(new = TRUE)
  plot(soft$Power, soft$mean.k., type = "l", axes = FALSE, xlab = "", ylab = "",
       xlim = c(1, 50), ylim = c(0, max(soft$mean.k., na.rm = TRUE) * 1.06),
       col = "#606060", lwd = 1.0, lty = 3)
  axis(4, col.axis = "#505050", cex.axis = 0.78)
  mtext("Mean connectivity", side = 4, line = 2.7, col = "#505050", cex = 0.82)
  panel_label("A")
}

draw_module_structure <- function() {
  old <- par(mar = c(4.5, 5.2, 3.0, 1.2), xaxs = "i", yaxs = "i")
  on.exit(par(old))
  all_modules <- c("turquoise", "blue", "brown", "grey")
  all_sizes <- c(turquoise = 4085, blue = 94, brown = 52, grey = 58)
  stopifnot(all(module_sizes$Freq[match(module_order, module_sizes$module_colors)] ==
                  all_sizes[module_order]))
  grey_size <- sum(network_object$network$colors == "grey")
  stopifnot(grey_size == all_sizes[["grey"]])
  values <- rev(unname(all_sizes[all_modules]))
  labels <- rev(all_modules)
  bars <- barplot(values, names.arg = labels, horiz = TRUE,
                  col = unname(module_palette[labels]), border = NA,
                  xlim = c(0, 4450), las = 1, cex.names = 0.92,
                  main = "Revised module size distribution",
                  xlab = "Number of genes", cex.main = 1.05,
                  axes = FALSE)
  axis(1, at = seq(0, 4000, by = 1000), cex.axis = 0.84)
  abline(v = seq(0, 4000, by = 1000), col = "#E8E8E8", lwd = 0.7)
  box(bty = "l", col = "#555555", lwd = 0.8)
  text(values + 55, bars, format(values, big.mark = ","), adj = c(0, 0.5),
       cex = 0.88, font = 2, col = "#252525", xpd = TRUE)
  panel_label("B")
}

draw_matrix <- function(mat, row_labels, col_labels, title, limit, labels = NULL,
                        row_colors = NULL, footnote = NULL) {
  nr <- nrow(mat); nc <- ncol(mat)
  plot(c(0.5, nc + 0.5), c(0.5, nr + 0.5), type = "n", axes = FALSE,
       xlab = "", ylab = "", main = title, cex.main = 0.98, xaxs = "i", yaxs = "i")
  for (i in seq_len(nr)) for (j in seq_len(nc)) {
    y <- nr - i + 1
    rect(j - 0.5, y - 0.5, j + 0.5, y + 0.5,
         col = diverging_color(mat[i, j], limit), border = "white", lwd = 1.2)
    txt <- if (is.null(labels)) sprintf("%.2f", mat[i, j]) else labels[i, j]
    text(j, y, txt, cex = if (is.null(labels)) 0.94 else 0.72)
  }
  axis(1, at = seq_len(nc), labels = col_labels, tick = FALSE, line = -0.6,
       cex.axis = 0.82)
  axis(2, at = nr:1, labels = row_labels, tick = FALSE, las = 1, line = -0.6,
       cex.axis = 0.86)
  if (!is.null(footnote)) mtext(footnote, side = 1, line = 2.2, cex = 0.61, adj = 0)
}

draw_eigengene_cor <- function() {
  old <- par(mar = c(5.0, 5.3, 3.0, 1.0))
  on.exit(par(old))
  me_names <- paste0("ME", module_order)
  mat <- cor(MEs[, me_names, drop = FALSE], use = "pairwise.complete.obs")
  draw_matrix(mat, module_order, module_order, "Module eigengene correlations", 1,
              row_colors = module_palette[module_order])
  panel_label("C")
}

draw_module_trait <- function() {
  old <- par(mar = c(6.0, 5.3, 3.0, 1.0))
  on.exit(par(old))
  mat <- matrix(NA_real_, nrow = length(module_order), ncol = length(effect_order),
                dimnames = list(module_order, effect_order))
  lab <- matrix("", nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  for (i in seq_len(nrow(mixed))) {
    mod <- sub("^ME", "", mixed$module[i])
    eff <- mixed$effect[i]
    if (mod %in% module_order && eff %in% effect_order) {
      mat[mod, eff] <- mixed$coefficient[i]
      lab[mod, eff] <- sprintf("b = %.2f\nP = %.3f\nFDR = %.3f",
                               mixed$coefficient[i], mixed$P_value[i], mixed$BH_FDR_global[i])
    }
  }
  draw_matrix(mat, module_order, effect_labels,
              "Subject-aware module-trait associations", max(abs(mat), na.rm = TRUE),
              labels = lab, row_colors = module_palette[module_order])
  panel_label("D")
}

draw_figure <- function() {
  layout(matrix(1:4, nrow = 2, byrow = TRUE), widths = c(1.05, 0.95), heights = c(1, 1))
  par(family = "sans", fg = "#222222", col.axis = "#222222", col.lab = "#222222",
      cex.axis = 0.90, cex.lab = 0.96, lend = "round")
  draw_soft_threshold()
  draw_module_structure()
  draw_eigengene_cor()
  draw_module_trait()
}

pdf_file <- file.path(output_dir, "Figure3_revised_WGCNA.pdf")
png_file <- file.path(output_dir, "Figure3_revised_WGCNA.png")
tif_file <- file.path(output_dir, "Figure3_revised_WGCNA.tiff")

cairo_pdf(pdf_file, width = 10.0, height = 8.2, family = "Arial")
draw_figure()
dev.off()

png(png_file, width = 10.0, height = 8.2, units = "in", res = 600,
    type = "cairo", family = "Arial", bg = "white")
draw_figure()
dev.off()

tiff(tif_file, width = 10.0, height = 8.2, units = "in", res = 600,
     compression = "lzw", type = "cairo", family = "Arial", bg = "white")
draw_figure()
dev.off()

legend_text <- paste(
  "Figure 3. Revised exploratory weighted gene co-expression network analysis of GSE75665.",
  "(A) Soft-threshold diagnostics. Powers 1-30 formed the prespecified candidate range; powers 31-50 were assessed only as a diagnostic extension.",
  "No power reached the prespecified scale-free topology fit threshold of R2 = 0.85. Beta = 30 was selected by the documented fallback rule (R2 = 0.103; mean connectivity = 753.02), and the criterion was not reached.",
  "(B) Module size distribution for the revised signed network: turquoise, 4,085 genes; blue, 94 genes; brown, 52 genes; and grey/unclassified, 58 genes.",
  "(C) Pearson correlations among the three biological module eigengenes.",
  "(D) Subject-aware mixed-effects associations of module eigengenes with AMS status, altitude status, and the AMS-by-altitude interaction.",
  "Each cell reports the model coefficient (b), nominal P value (P), and global Benjamini-Hochberg-adjusted P value (FDR) across all nine module-effect tests.",
  "No module-effect association remained significant at BH-adjusted P < 0.05."
)
writeLines(legend_text, file.path(output_dir, "Figure3_revised_WGCNA_legend.txt"))

new_sheets <- c("S4A_WGCNA_parameters", "S4B_gene_modules", "S4C_module_sizes",
                "S4D_module_trait", "S4E_gene_metrics", "S4F_candidates",
                "S4G_soft_threshold", "S4H_sample_QC")

config <- read_result("WGCNA_analysis_config.csv")
filter_summary <- read_result("WGCNA_gene_filter_summary.csv")
parameters <- data.frame(
  category = c(rep("Analysis configuration", nrow(config)), rep("Gene filtering", nrow(filter_summary))),
  parameter = c(config$item, filter_summary$stage),
  value = c(config$value, filter_summary$genes),
  detail = c(rep("", nrow(config)), filter_summary$rule),
  source_file = c(rep("WGCNA_analysis_config.csv", nrow(config)),
                  rep("WGCNA_gene_filter_summary.csv", nrow(filter_summary)))
)

tables <- list(
  S4A_WGCNA_parameters = parameters,
  S4B_gene_modules = read_result("WGCNA_module_assignments.csv"),
  S4C_module_sizes = read_result("WGCNA_module_sizes.csv"),
  S4D_module_trait = mixed,
  S4E_gene_metrics = read_result("WGCNA_gene_module_metrics.csv"),
  S4F_candidates = read_result("WGCNA_candidate_gene_metrics.csv"),
  S4G_soft_threshold = soft,
  S4H_sample_QC = read_result("WGCNA_sample_QC.csv")
)

# CSV mirrors make each supplementary result independently inspectable.
for (s in names(tables)) {
  write.csv(tables[[s]], file.path(output_dir, paste0(s, ".csv")), row.names = FALSE)
}

cat("Figure PDF:", pdf_file, "\n")
cat("Figure PNG:", png_file, "\n")
cat("Figure TIFF:", tif_file, "\n")
