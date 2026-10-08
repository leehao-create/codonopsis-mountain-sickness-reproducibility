#!/usr/bin/env Rscript

# Format archived manuscript-primary WGCNA outputs into release tables.
# This script does not rebuild the network or alter module assignments.

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
primary <- file.path(root, "analysis", "WGCNA", "manuscript_primary")
reviewer <- file.path(root, "analysis", "WGCNA", "reviewer_checks")
reference <- file.path(root, "analysis", "WGCNA", "data", "reference", "Codonopsis_54_overlap_targets.csv")

label_colors <- c(
  `0` = "grey", `1` = "turquoise", `2` = "blue", `3` = "brown",
  `4` = "yellow", `5` = "green", `6` = "red", `7` = "black",
  `8` = "pink", `9` = "magenta", `10` = "purple", `11` = "greenyellow"
)

numeric <- read.csv(file.path(primary, "results", "WGCNA_module_numeric_labels.csv"),
                    row.names = 1, check.names = FALSE)
assignments <- data.frame(
  gene = rownames(numeric),
  module_label = numeric[[1]],
  module_color = unname(label_colors[as.character(numeric[[1]])]),
  stringsAsFactors = FALSE
)
stopifnot(nrow(assignments) == 2810L, !anyNA(assignments$module_color))
write.csv(assignments, file.path(primary, "results", "WGCNA_module_assignments.csv"), row.names = FALSE)

sizes <- as.data.frame(table(assignments$module_label, assignments$module_color), stringsAsFactors = FALSE)
sizes <- sizes[sizes$Freq > 0, ]
names(sizes) <- c("module_label", "module_color", "module_size")
sizes <- sizes[order(as.integer(sizes$module_label)), ]
write.csv(sizes, file.path(primary, "results", "WGCNA_module_sizes.csv"), row.names = FALSE)

target_table <- read.csv(reference, check.names = FALSE)
target_col <- names(target_table)[[1]]
target_map <- merge(data.frame(gene = unique(target_table[[target_col]])), assignments,
                    by = "gene", all.x = TRUE, sort = FALSE)
target_map$in_retained_gene_set <- !is.na(target_map$module_label)
write.csv(target_map, file.path(reviewer, "results", "WGCNA_candidate_target_module_mapping.csv"), row.names = FALSE)

p4 <- read.csv(file.path(primary, "results", "WGCNA_module_trait_4group_pvalues.csv"),
               row.names = 1, check.names = FALSE)
r4 <- read.csv(file.path(primary, "results", "WGCNA_module_trait_4group_correlations.csv"),
               row.names = 1, check.names = FALSE)
stopifnot(identical(dim(p4), c(12L, 4L)), identical(dim(r4), c(12L, 4L)))
fdr4 <- matrix(p.adjust(as.vector(as.matrix(p4)), method = "BH"),
               nrow = nrow(p4), dimnames = dimnames(p4))
long4 <- do.call(rbind, lapply(seq_len(nrow(p4)), function(i) {
  data.frame(module = rownames(p4)[i], trait = colnames(p4),
             Pearson_r = as.numeric(r4[i, ]), nominal_P = as.numeric(p4[i, ]),
             BH_adjusted_P_global_48 = as.numeric(fdr4[i, ]), check.names = FALSE)
}))
write.csv(long4, file.path(reviewer, "results", "WGCNA_4group_48test_BH_results.csv"), row.names = FALSE)

binary <- read.csv(file.path(reviewer, "results", "WGCNA_module_AMS_statistics.csv"), check.names = FALSE)
binary$BH_adjusted_P_12_tests <- p.adjust(binary$P_value, method = "BH")
write.csv(binary, file.path(reviewer, "results", "WGCNA_binary_AMS_statistics_with_BH.csv"), row.names = FALSE)

cat("Formatted", nrow(assignments), "module assignments; minimum 48-test BH =",
    min(long4$BH_adjusted_P_global_48), "; minimum binary-AMS BH =",
    min(binary$BH_adjusted_P_12_tests), "\n")
