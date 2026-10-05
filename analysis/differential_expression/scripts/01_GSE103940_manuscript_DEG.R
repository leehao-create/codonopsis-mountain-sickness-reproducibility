#!/usr/bin/env Rscript

# Reproduce the GSE103940 differential-expression table reported in the
# submitted manuscript. This script preserves the archived analysis exactly;
# it does not substitute the later count-based sensitivity analysis.

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
  getwd()
}
module_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

suppressPackageStartupMessages(library(limma))

input_file <- file.path(
  module_dir, "data", "processed", "GSE103940_FPKM_gene_symbol_matrix.csv"
)
metadata_file <- file.path(
  module_dir, "data", "metadata", "GSE103940_subject_pairing.csv"
)
output_dir <- Sys.getenv(
  "DE_OUTPUT_DIR",
  file.path(module_dir, "reproduced_results", "GSE103940_manuscript")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

expression <- read.csv(input_file, row.names = 1, check.names = FALSE)
metadata <- read.csv(metadata_file, check.names = FALSE)

# The archived matrix columns retain the GEO supplementary-file names.
stopifnot(
  nrow(metadata) == 22L,
  all(table(metadata$subject_id) == 2L),
  identical(colnames(expression), metadata$source_file)
)

# Frozen manuscript workflow: remove genes with zero total expression, apply
# log2 to FPKM values, and replace log2(0) with zero. No additional
# normalization was recorded in the archived workflow.
expression <- expression[rowSums(expression) != 0, , drop = FALSE]
expression <- log2(as.matrix(expression))
expression[is.infinite(expression) & expression < 0] <- 0

# This is the design used by the archived manuscript analysis. Subject was not
# included in this model. The stored coefficient is plain minus high altitude;
# manuscript directions are therefore the inverse of the logFC sign.
condition <- factor(metadata$condition, levels = c("high_altitude", "plain"))
design <- model.matrix(~ 0 + condition)
colnames(design) <- c("high", "plain")
rownames(design) <- metadata$source_file
contrast <- makeContrasts(plain - high, levels = design)

fit <- lmFit(expression, design)
fit <- contrasts.fit(fit, contrast)
fit <- eBayes(fit)
deg <- topTable(fit, coef = 1, number = Inf, sort.by = "logFC")
deg <- na.omit(deg)

primary <- deg$adj.P.Val < 0.05 & abs(deg$logFC) > 1.5
high_altitude_up <- primary & deg$logFC < -1.5
high_altitude_down <- primary & deg$logFC > 1.5
stopifnot(
  sum(primary) == 2782L,
  sum(high_altitude_up) == 118L,
  sum(high_altitude_down) == 2664L
)

write.csv(
  deg,
  file.path(output_dir, "GSE103940_DEG_all_plain_minus_high.csv")
)
write.csv(
  deg[primary, , drop = FALSE],
  file.path(output_dir, "GSE103940_DEG_significant_plain_minus_high.csv")
)
write.csv(
  data.frame(
    reported_comparison = "high altitude - plain",
    adjusted_P_threshold = 0.05,
    absolute_log2FC_threshold = 1.5,
    total_DEGs = sum(primary),
    up_in_high_altitude = sum(high_altitude_up),
    down_in_high_altitude = sum(high_altitude_down),
    stored_table_contrast = "plain - high altitude"
  ),
  file.path(output_dir, "GSE103940_manuscript_DEG_counts.csv"),
  row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

