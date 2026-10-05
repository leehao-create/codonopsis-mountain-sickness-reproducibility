options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
module_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
suppressPackageStartupMessages({library(Matrix); library(ggplot2)})

sc_out <- Sys.getenv("SINGLE_CELL_OUTPUT_DIR", file.path(module_dir, "reproduced_results"))
dir.create(sc_out, recursive = TRUE, showWarnings = FALSE)
target_file <- file.path(module_dir, "..", "network_pharmacology", "data", "derived",
                         "Codonopsis_mountain_sickness_54_targets.csv")
targets <- unique(trimws(read.csv(target_file, check.names = FALSE)[[1]]))
stopifnot(length(targets) == 54L, !anyDuplicated(targets))
plot_genes <- c("ATM", "SNCA", "BCL2", "APP", "HBB", "CA2", "GATA1")

external_dir <- Sys.getenv("HCL_RDS_DIR", file.path(module_dir, "data", "external"))
objects <- c(blood = file.path(external_dir, "blood1.rds"),
             lung = file.path(external_dir, "lung1.rds"))
if (!all(file.exists(objects))) {
  stop("The Human Cell Landscape RDS objects are not redistributed. Place blood1.rds ",
       "and lung1.rds in HCL_RDS_DIR; see README.md and external_data_manifest.csv.")
}

expression_rows <- list()
coverage_rows <- list()
dataset_rows <- list()
donor_rows <- list()
counter <- 0L
for (tissue in names(objects)) {
  object <- readRDS(objects[[tissue]])
  metadata <- attr(object, "meta.data")
  assay <- attr(object, "assays")[["RNA"]]
  counts <- attr(assay, "counts")
  normalized <- attr(assay, "data")
  stopifnot(ncol(counts) == nrow(metadata), ncol(normalized) == nrow(metadata),
            identical(colnames(counts), rownames(metadata)),
            "celltype" %in% names(metadata), "donor_id" %in% names(metadata),
            identical(rownames(counts), rownames(normalized)))
  metadata$celltype <- as.character(metadata$celltype)
  present <- targets %in% rownames(normalized)
  coverage_rows[[tissue]] <- data.frame(
    tissue = tissue, gene = targets, gene_present = present,
    any_expression = vapply(targets, function(gene) {
      gene %in% rownames(counts) && sum(counts[gene, ]) > 0
    }, logical(1)), stringsAsFactors = FALSE
  )
  for (celltype in sort(unique(metadata$celltype))) {
    idx <- which(metadata$celltype == celltype)
    genes <- targets[present]
    avg <- Matrix::rowMeans(normalized[genes, idx, drop = FALSE])
    pct <- Matrix::rowMeans(counts[genes, idx, drop = FALSE] > 0) * 100
    counter <- counter + 1L
    expression_rows[[counter]] <- data.frame(
      tissue = tissue, celltype = celltype, n_cells = length(idx), gene = genes,
      average_log_normalized_expression = as.numeric(avg),
      percent_expressed = as.numeric(pct), stringsAsFactors = FALSE
    )
  }
  donor_table <- as.data.frame(table(donor_id = as.character(metadata$donor_id)),
                               stringsAsFactors = FALSE)
  donor_table <- donor_table[donor_table$Freq > 0L, ]
  donor_table$tissue <- tissue
  donor_rows[[tissue]] <- donor_table[, c("tissue", "donor_id", "Freq")]
  names(donor_rows[[tissue]])[3] <- "cell_count"
  tissue_stage <- if ("tissue_original" %in% names(metadata)) {
    paste(names(sort(table(as.character(metadata$tissue_original)), decreasing = TRUE)),
          as.integer(sort(table(as.character(metadata$tissue_original)), decreasing = TRUE)),
          sep = ":", collapse = "; ")
  } else "Unavailable"
  dataset_rows[[tissue]] <- data.frame(
    tissue = tissue, cells = nrow(metadata), genes = nrow(normalized),
    donors = length(unique(as.character(metadata$donor_id))),
    batches = length(unique(as.character(metadata$batch))),
    celltypes = length(unique(metadata$celltype)),
    developmental_composition = tissue_stage,
    disease_metadata = paste(unique(as.character(metadata$disease)), collapse = ";"),
    assay_metadata = paste(unique(as.character(metadata$assay)), collapse = ";"),
    object_file = objects[[tissue]], stringsAsFactors = FALSE
  )
}

expression <- do.call(rbind, expression_rows)
coverage <- do.call(rbind, coverage_rows)
dataset_summary <- do.call(rbind, dataset_rows)
donor_summary <- do.call(rbind, donor_rows)

major_rows <- list(); k <- 0L
for (tissue in names(objects)) {
  for (gene in targets) {
    z <- expression[expression$tissue == tissue & expression$gene == gene, ]
    k <- k + 1L
    if (!nrow(z)) {
      major_rows[[k]] <- data.frame(tissue = tissue, gene = gene,
        gene_present = FALSE, any_expression = FALSE, major_celltype = NA_character_,
        major_celltype_percent_expressed = NA_real_,
        major_celltype_average_expression = NA_real_)
    } else {
      z <- z[order(-z$percent_expressed, -z$average_log_normalized_expression,
                   z$celltype), ]
      major_rows[[k]] <- data.frame(tissue = tissue, gene = gene,
        gene_present = TRUE, any_expression = any(z$percent_expressed > 0),
        major_celltype = z$celltype[1],
        major_celltype_percent_expressed = z$percent_expressed[1],
        major_celltype_average_expression = z$average_log_normalized_expression[1])
    }
  }
}
major <- do.call(rbind, major_rows)

write.csv(expression, file.path(sc_out, "single_cell_all54_expression_by_celltype.csv"),
          row.names = FALSE)
write.csv(expression[expression$gene %in% plot_genes, ],
          file.path(sc_out, "single_cell_revised_candidate_expression_by_celltype.csv"),
          row.names = FALSE)
write.csv(coverage, file.path(sc_out, "single_cell_all54_gene_coverage.csv"),
          row.names = FALSE)
write.csv(major, file.path(sc_out, "single_cell_all54_major_celltypes.csv"),
          row.names = FALSE)
write.csv(dataset_summary, file.path(sc_out, "single_cell_dataset_summary.csv"),
          row.names = FALSE)
write.csv(donor_summary, file.path(sc_out, "single_cell_donor_cell_counts.csv"),
          row.names = FALSE)

plot_data <- expression[expression$gene %in% plot_genes, ]
plot_data$gene <- factor(plot_data$gene, levels = plot_genes)
plot_data$tissue <- factor(plot_data$tissue, levels = c("blood", "lung"),
                           labels = c("Reference blood", "Reference lung"))
cell_order <- c(sort(unique(plot_data$celltype[plot_data$tissue == "Reference blood"])),
                sort(unique(plot_data$celltype[plot_data$tissue == "Reference lung"])))
plot_data$celltype <- factor(plot_data$celltype, levels = unique(cell_order))

p <- ggplot(plot_data, aes(gene, celltype)) +
  geom_point(aes(size = percent_expressed,
                 colour = average_log_normalized_expression), alpha = .92) +
  facet_grid(tissue ~ ., scales = "free_y", space = "free_y") +
  scale_size_continuous(range = c(.2, 5.8), breaks = c(0, 25, 50, 75, 100),
                        limits = c(0, 100), name = "Percent\nexpressed") +
  scale_colour_gradientn(colours = c("#F2F2F2", "#74A9CF", "#045A8D"),
                         name = "Average log-normalized\nexpression") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, colour = "black"),
        axis.text.y = element_text(colour = "black"),
        strip.background = element_rect(fill = "#F2F2F2", colour = "#666666",
                                        linewidth = .4),
        strip.text = element_text(face = "bold"),
        axis.line = element_blank(), axis.ticks = element_blank(),
        panel.border = element_rect(fill = NA, colour = "#777777", linewidth = .35),
        panel.spacing = grid::unit(.10, "inches"), legend.position = "right",
        plot.margin = margin(6, 8, 6, 6))

ggsave(file.path(sc_out, "candidate_gene_celltype_dotplot.pdf"), p,
       width = 7.2, height = 7.0, units = "in", device = cairo_pdf,
       family = "Arial")
ggsave(file.path(sc_out, "candidate_gene_celltype_dotplot.png"), p,
       width = 7.2, height = 7.0, units = "in", dpi = 600, type = "cairo-png")
system2("mogrify", c("-units", "PixelsPerInch", "-density", "600",
                      shQuote(file.path(sc_out, "candidate_gene_celltype_dotplot.png"))))

input_files <- c(unname(objects), target_file)
write.csv(data.frame(path = input_files, md5 = unname(tools::md5sum(input_files))),
          file.path(sc_out, "single_cell_input_md5.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(sc_out, "single_cell_sessionInfo.txt"))

print(dataset_summary)
print(major[major$gene %in% plot_genes, ])
