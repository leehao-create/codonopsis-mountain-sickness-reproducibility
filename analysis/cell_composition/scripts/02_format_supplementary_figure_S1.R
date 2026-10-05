#!/usr/bin/env Rscript

# Layout-only regeneration of Supplementary Figure S1.
# This script reads frozen marker-score results and does not recompute statistics.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
} else {
  getwd()
}

module_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
input_file <- file.path(module_dir, "results", "cell_type_marker_score_results.csv")
output_dir <- Sys.getenv("CELL_COMPOSITION_FIGURE_DIR", file.path(module_dir, "reproduced_figures"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(input_file))
results <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)

comparison_order <- c(
  "paired high altitude - plain marker-score change",
  "paired high altitude - plain marker-score change; non-globin-factor TMM sensitivity",
  "average paired high altitude - plain marker-score change",
  "AMS versus non-AMS difference in paired marker-score change"
)
panel_titles <- c(
  "High altitude \u2212 plain marker-score change",
  "High altitude \u2212 plain marker-score change (non-globin sensitivity)",
  "Mean paired marker-score change",
  "AMS vs non-AMS difference in paired marker-score change"
)
panel_labels <- LETTERS[1:4]
lineage_order <- c(
  "T_cells", "B_cells", "NK_cells", "Neutrophils", "Monocytes",
  "Erythroid_reticulocyte"
)

plot_df <- results[
  results$comparison %in% comparison_order & results$dataset != "GSE260910",
]
stopifnot(
  nrow(plot_df) == length(comparison_order) * length(lineage_order),
  all(is.finite(plot_df$effect)),
  all(is.finite(plot_df$CI_lower)),
  all(is.finite(plot_df$CI_upper)),
  all(is.finite(plot_df$BH_FDR))
)

significant_colour <- "#B33A3A"
neutral_colour <- "#52606D"
zero_colour <- "#C7CCD1"
axis_colour <- "#262B30"

plot_figure <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)

  par(
    mfrow = c(4, 1),
    mar = c(6.1, 11.2, 4.4, 2.0),
    oma = c(4.4, 0.5, 0.5, 0.3),
    family = "Arial",
    ps = 10,
    mgp = c(3.5, 0.75, 0),
    tcl = -0.25,
    las = 1,
    xaxs = "r",
    yaxs = "i",
    lend = "round",
    bg = "white",
    fg = axis_colour,
    col.axis = axis_colour,
    col.lab = axis_colour,
    col.main = axis_colour
  )

  for (i in seq_along(comparison_order)) {
    z <- plot_df[plot_df$comparison == comparison_order[[i]], ]
    z <- z[match(lineage_order, z$cell_type), ]
    stopifnot(!anyNA(z$cell_type))

    # Reverse y coordinates so the first lineage appears at the top.
    y <- rev(seq_len(nrow(z)))
    limits <- range(c(z$CI_lower, z$CI_upper, 0), finite = TRUE)
    pad <- diff(limits) * 0.08
    if (!is.finite(pad) || pad == 0) pad <- 0.1
    x_limits <- limits + c(-pad, pad)
    point_colours <- ifelse(z$BH_FDR < 0.05, significant_colour, neutral_colour)

    plot(
      NA_real_, NA_real_,
      xlim = x_limits,
      ylim = c(0.5, nrow(z) + 0.5),
      xaxt = "n", yaxt = "n",
      xlab = "", ylab = "",
      bty = "l",
      lwd = 0.9
    )
    abline(v = 0, col = zero_colour, lwd = 1.0)
    segments(
      x0 = z$CI_lower, y0 = y,
      x1 = z$CI_upper, y1 = y,
      col = point_colours, lwd = 1.7
    )
    points(z$effect, y, pch = 19, cex = 1.25, col = point_colours)

    axis(1, cex.axis = 0.95, lwd = 0.8, lwd.ticks = 0.8)
    axis(
      2, at = y, labels = z$cell_type,
      cex.axis = 0.95, lwd = 0.8, lwd.ticks = 0.8,
      las = 1, padj = 0.15
    )
    title(
      main = panel_titles[[i]],
      xlab = "Mean change in aggregate marker z-score (95% CI)",
      line = 3.6,
      cex.main = 1.15,
      cex.lab = 1.05,
      font.main = 2,
      font.lab = 1
    )
    mtext(
      panel_labels[[i]], side = 3, line = 2.65,
      adj = -0.085, cex = 1.2, font = 2, xpd = NA
    )
  }

  mtext(
    "Expression-based lineage marker proxy; not an estimated cell proportion",
    side = 1, outer = TRUE, line = 2.35, cex = 0.95, font = 2,
    col = axis_colour
  )
}

pdf_file <- file.path(
  output_dir, "Supplementary_Figure_S1_Cell_Composition_Marker_Audit.pdf"
)
png_file <- file.path(
  output_dir, "Supplementary_Figure_S1_Cell_Composition_Marker_Audit.png"
)

# 180 mm manuscript width; extra height is intentional to preserve readable type.
cairo_pdf(
  filename = pdf_file,
  width = 180 / 25.4,
  height = 390 / 25.4,
  family = "Arial",
  onefile = TRUE,
  bg = "white"
)
plot_figure()
dev.off()

png(
  filename = png_file,
  width = 180 / 25.4,
  height = 390 / 25.4,
  units = "in",
  res = 600,
  type = "cairo",
  family = "Arial",
  bg = "white"
)
plot_figure()
dev.off()

message("Input:  ", input_file)
message("PDF:    ", pdf_file)
message("PNG:    ", png_file)
