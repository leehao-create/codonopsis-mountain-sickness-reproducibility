options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
module_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

out <- Sys.getenv("CELL_COMPOSITION_OUTPUT_DIR", file.path(module_dir, "reproduced_results"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Compact canonical lineage-marker panel. These genes are used as expression
# proxies only; their abundance cannot identify cell proportions uniquely.
marker_sets <- list(
  T_cells = c("CD3D", "CD3E", "CD3G", "TRAC", "LCK", "IL7R"),
  B_cells = c("CD19", "MS4A1", "CD79A", "CD79B", "CD37", "CD74"),
  NK_cells = c("NKG7", "GNLY", "KLRD1", "PRF1", "GZMB", "NCAM1"),
  Neutrophils = c("FCGR3B", "CSF3R", "CEACAM8", "S100A8", "S100A9", "FPR1"),
  Monocytes = c("LST1", "LILRB1", "CTSS", "FCN1", "TYMP", "CTSD"),
  Erythroid_reticulocyte = c("HBB", "HBA1", "HBA2", "ALAS2", "AHSP",
                             "GYPA", "SLC4A1", "CA1", "EPB42", "BPGM")
)
marker_panel <- do.call(rbind, lapply(names(marker_sets), function(cell_type) {
  data.frame(cell_type = cell_type, gene = marker_sets[[cell_type]],
             panel_role = "canonical lineage-expression proxy; not a unique proportion marker")
}))
write.csv(marker_panel, file.path(out, "cell_type_marker_panel.csv"), row.names = FALSE)

direction_label <- function(effect, positive, negative) {
  ifelse(is.na(effect), "not_estimable",
         ifelse(effect > 0, positive, ifelse(effect < 0, negative, "no_change")))
}

make_scores <- function(expression, sets) {
  scores <- matrix(NA_real_, nrow = length(sets), ncol = ncol(expression),
                   dimnames = list(names(sets), colnames(expression)))
  coverage <- data.frame()
  for (cell_type in names(sets)) {
    present <- intersect(sets[[cell_type]], rownames(expression))
    z <- t(scale(t(expression[present, , drop = FALSE])))
    z[!is.finite(z)] <- NA_real_
    scores[cell_type, ] <- colMeans(z, na.rm = TRUE)
    coverage <- rbind(coverage, data.frame(cell_type = cell_type,
                                           requested_markers = length(sets[[cell_type]]),
                                           present_markers = length(present),
                                           genes_present = paste(present, collapse = ";"),
                                           genes_absent = paste(setdiff(sets[[cell_type]], present),
                                                                collapse = ";")))
  }
  list(scores = scores, coverage = coverage)
}

paired_deltas <- function(scores, metadata, subject_col, condition_col,
                          plain = "plain", high = "high_altitude") {
  subjects <- unique(metadata[[subject_col]])
  delta <- sapply(subjects, function(subject) {
    idx <- which(metadata[[subject_col]] == subject)
    hi <- idx[metadata[[condition_col]][idx] == high]
    lo <- idx[metadata[[condition_col]][idx] == plain]
    stopifnot(length(hi) == 1L, length(lo) == 1L)
    scores[, hi] - scores[, lo]
  })
  colnames(delta) <- subjects
  delta
}

one_sample_summary <- function(delta, dataset, comparison, direction_positive,
                               direction_negative) {
  ans <- do.call(rbind, lapply(rownames(delta), function(cell_type) {
    x <- as.numeric(delta[cell_type, ])
    test <- t.test(x, mu = 0)
    data.frame(dataset = dataset, comparison = comparison, cell_type = cell_type,
               marker_count = NA_integer_, effect = mean(x),
               CI_lower = unname(test$conf.int[1]), CI_upper = unname(test$conf.int[2]),
               P_value = test$p.value)
  }))
  ans$BH_FDR <- p.adjust(ans$P_value, method = "BH")
  ans$direction <- direction_label(ans$effect, direction_positive, direction_negative)
  ans$significance <- ifelse(ans$BH_FDR < 0.05, "BH-FDR < 0.05", "not significant")
  ans$composition_interpretation <- ifelse(
    ans$BH_FDR < 0.05,
    "coherent lineage-marker shift; consistent with, but does not prove, a composition change",
    "no statistically supported shift in the aggregate marker proxy")
  ans
}

# GSE103940: frozen paired count/TMM/voom object used only for the independent
# reviewer-requested lineage-marker audit. It is not the manuscript-locked
# GSE103940 primary DEG table (2,782 genes from the archived FPKM workflow).
fit103_path <- file.path(module_dir, "data", "dependencies", "GSE103940_final_A_B_fit.rds")
deg103_path <- file.path(module_dir, "data", "dependencies",
                         "GSE103940_cell_composition_count_voom_DEG_all.csv")
deg103_b_path <- file.path(module_dir, "data", "dependencies", "GSE103940_analysis_B_DEG_all.csv")
fit103 <- readRDS(fit103_path)
deg103 <- read.csv(deg103_path, check.names = FALSE)
deg103_b <- read.csv(deg103_b_path, check.names = FALSE)
expr103 <- fit103$voom_A$E
meta103 <- fit103$metadata
stopifnot(identical(colnames(expr103), meta103$gsm), nrow(meta103) == 22L,
          all(table(meta103$subject_id) == 2L))

gene103 <- merge(marker_panel, deg103[, c("gene", "logFC", "P.Value", "adj.P.Val")],
                 by = "gene", all.x = TRUE, sort = FALSE)
gene103 <- gene103[match(marker_panel$gene, gene103$gene), ]
gene103$dataset <- "GSE103940"
gene103$comparison <- "high altitude - plain (paired TMM/voom-limma)"
gene103$input_available <- gene103$gene %in% rownames(fit103$y_A$counts)
gene103$expression_eligible <- gene103$gene %in% rownames(expr103)
gene103$direction <- direction_label(gene103$logFC, "higher_in_high_altitude",
                                     "lower_in_high_altitude")
gene103$significance <- ifelse(is.na(gene103$adj.P.Val), "not expression-eligible",
                               ifelse(gene103$adj.P.Val < 0.05, "BH-FDR < 0.05",
                                      "not significant"))
gene103$interpretation <- "individual marker expression; compatible with both abundance and cell-intrinsic effects"
gene103$whether_supports_composition_shift <- ifelse(
  gene103$adj.P.Val < 0.05,
  "consistent direction for this marker only; aggregate panel required",
  "no individual-marker statistical support")
gene103 <- gene103[, c("dataset", "comparison", "cell_type", "gene", "input_available",
                       "expression_eligible", "logFC", "P.Value", "adj.P.Val",
                       "direction", "significance", "interpretation",
                       "whether_supports_composition_shift")]

gene103_b <- merge(marker_panel, deg103_b[, c("gene", "logFC", "P.Value", "adj.P.Val")],
                   by = "gene", all.x = TRUE, sort = FALSE)
gene103_b <- gene103_b[match(marker_panel$gene, gene103_b$gene), ]
gene103_b$dataset <- "GSE103940"
gene103_b$comparison <- "high altitude - plain (non-globin-factor TMM sensitivity)"
gene103_b$input_available <- gene103_b$gene %in% rownames(fit103$y_B$counts)
gene103_b$expression_eligible <- gene103_b$gene %in% rownames(fit103$voom_B$E)
gene103_b$direction <- direction_label(gene103_b$logFC, "higher_in_high_altitude",
                                       "lower_in_high_altitude")
gene103_b$significance <- ifelse(is.na(gene103_b$adj.P.Val), "not expression-eligible",
                                 ifelse(gene103_b$adj.P.Val < 0.05, "BH-FDR < 0.05",
                                        "not significant"))
gene103_b$interpretation <- "globin-composition sensitivity model; individual marker expression"
gene103_b$whether_supports_composition_shift <- ifelse(
  gene103_b$adj.P.Val < 0.05,
  "consistent direction for this marker only; aggregate panel required",
  "no individual-marker statistical support")
gene103_b <- gene103_b[, colnames(gene103)]

scores103 <- make_scores(expr103, marker_sets)
delta103 <- paired_deltas(scores103$scores, meta103, "subject_id", "condition")
score103 <- one_sample_summary(delta103, "GSE103940",
                               "paired high altitude - plain marker-score change",
                               "higher_in_high_altitude", "lower_in_high_altitude")
score103$marker_count <- scores103$coverage$present_markers[match(score103$cell_type,
                                                                  scores103$coverage$cell_type)]

scores103_b <- make_scores(fit103$voom_B$E, marker_sets)
delta103_b <- paired_deltas(scores103_b$scores, meta103, "subject_id", "condition")
score103_b <- one_sample_summary(
  delta103_b, "GSE103940",
  "paired high altitude - plain marker-score change; non-globin-factor TMM sensitivity",
  "higher_in_high_altitude", "lower_in_high_altitude")
score103_b$marker_count <- scores103_b$coverage$present_markers[
  match(score103_b$cell_type, scores103_b$coverage$cell_type)]

# GSE75665: frozen repeated-measures count/TMM/voom object.
fit756_path <- file.path(module_dir, "data", "dependencies",
                         "GSE75665_NCBI_counts_AMS_interaction_fit.rds")
interaction756_path <- file.path(module_dir, "data", "dependencies",
                                 "GSE75665_NCBI_counts_AMS_interaction_all_genes.csv")
fit756 <- readRDS(fit756_path)
interaction756 <- read.csv(interaction756_path, check.names = FALSE)
expr756 <- fit756$voom$E
meta756 <- fit756$metadata
delta756 <- fit756$subject_delta
stopifnot(identical(colnames(expr756), meta756$gsm), nrow(meta756) == 20L,
          all(table(meta756$subject_id) == 2L), ncol(delta756) == 10L)

# A transparent repeated-measures audit of the average altitude change: direct
# subject-level high-minus-plain deltas, tested against zero for all genes.
avg_effect <- rowMeans(delta756)
avg_sd <- apply(delta756, 1, sd)
avg_t <- avg_effect / (avg_sd / sqrt(ncol(delta756)))
avg_p <- 2 * pt(-abs(avg_t), df = ncol(delta756) - 1L)
avg_p[!is.finite(avg_p)] <- 1
avg_fdr <- p.adjust(avg_p, method = "BH")
avg756 <- data.frame(gene = rownames(delta756), logFC = avg_effect,
                     P.Value = avg_p, adj.P.Val = avg_fdr)

gene756_avg <- merge(marker_panel, avg756, by = "gene", all.x = TRUE, sort = FALSE)
gene756_avg <- gene756_avg[match(marker_panel$gene, gene756_avg$gene), ]
gene756_avg$dataset <- "GSE75665"
gene756_avg$comparison <- "average high altitude - plain change across 10 paired subjects"
gene756_avg$input_available <- gene756_avg$gene %in% rownames(expr756)
gene756_avg$expression_eligible <- gene756_avg$gene %in% rownames(delta756)
gene756_avg$direction <- direction_label(gene756_avg$logFC, "higher_in_high_altitude",
                                         "lower_in_high_altitude")
gene756_avg$significance <- ifelse(is.na(gene756_avg$adj.P.Val), "not expression-eligible",
                                   ifelse(gene756_avg$adj.P.Val < 0.05, "BH-FDR < 0.05",
                                          "not significant"))
gene756_avg$interpretation <- "direct mean of subject-level high-minus-plain voom changes; audit contrast"
gene756_avg$whether_supports_composition_shift <- ifelse(
  gene756_avg$adj.P.Val < 0.05,
  "consistent direction for this marker only; aggregate panel required",
  "no individual-marker statistical support")
gene756_avg <- gene756_avg[, colnames(gene103)]

gene756_int <- merge(marker_panel,
                     interaction756[, c("gene", "logFC", "P.Value", "adj.P.Val",
                                        "expression_eligible")],
                     by = "gene", all.x = TRUE, sort = FALSE)
gene756_int <- gene756_int[match(marker_panel$gene, gene756_int$gene), ]
gene756_int$dataset <- "GSE75665"
gene756_int$comparison <- "(AMS high - plain) - (non-AMS high - plain)"
gene756_int$input_available <- gene756_int$gene %in% interaction756$gene
gene756_int$direction <- direction_label(gene756_int$logFC,
                                         "larger_altitude_change_in_AMS",
                                         "smaller_altitude_change_in_AMS")
gene756_int$significance <- ifelse(is.na(gene756_int$adj.P.Val), "not expression-eligible",
                                   ifelse(gene756_int$adj.P.Val < 0.05, "BH-FDR < 0.05",
                                          "not significant"))
gene756_int$interpretation <- "frozen subject-aware AMS-by-altitude interaction"
gene756_int$whether_supports_composition_shift <- ifelse(
  gene756_int$adj.P.Val < 0.05,
  "consistent with an AMS-specific marker shift, but not proof of proportion change",
  "no AMS-specific individual-marker statistical support")
gene756_int <- gene756_int[, colnames(gene103)]

scores756 <- make_scores(expr756, marker_sets)
delta_score756 <- paired_deltas(scores756$scores, meta756, "subject_id", "altitude_status")
score756_avg <- one_sample_summary(delta_score756, "GSE75665",
                                   "average paired high altitude - plain marker-score change",
                                   "higher_in_high_altitude", "lower_in_high_altitude")
score756_avg$marker_count <- scores756$coverage$present_markers[
  match(score756_avg$cell_type, scores756$coverage$cell_type)]

subject_status <- vapply(colnames(delta_score756), function(subject) {
  unique(meta756$ams_status[meta756$subject_id == subject])
}, character(1))
score756_int <- do.call(rbind, lapply(rownames(delta_score756), function(cell_type) {
  x <- delta_score756[cell_type, subject_status == "AMS"]
  y <- delta_score756[cell_type, subject_status == "non_AMS"]
  test <- t.test(x, y)
  data.frame(dataset = "GSE75665",
             comparison = "AMS versus non-AMS difference in paired marker-score change",
             cell_type = cell_type,
             marker_count = scores756$coverage$present_markers[
               match(cell_type, scores756$coverage$cell_type)],
             effect = mean(x) - mean(y), CI_lower = unname(test$conf.int[1]),
             CI_upper = unname(test$conf.int[2]), P_value = test$p.value)
}))
score756_int$BH_FDR <- p.adjust(score756_int$P_value, method = "BH")
score756_int$direction <- direction_label(score756_int$effect,
                                          "larger_altitude_change_in_AMS",
                                          "smaller_altitude_change_in_AMS")
score756_int$significance <- ifelse(score756_int$BH_FDR < 0.05, "BH-FDR < 0.05",
                                    "not significant")
score756_int$composition_interpretation <- ifelse(
  score756_int$BH_FDR < 0.05,
  "coherent AMS-specific lineage-marker shift; consistent with, but does not prove, composition change",
  "no statistically supported AMS-specific shift in the aggregate marker proxy")

# GSE260910: descriptive values only. Disease and reported batch are identical,
# so no inferential P value, FDR, or cell-composition conclusion is estimable.
matrix260_path <- file.path(module_dir, "data", "dependencies", "GSE260910_processed_data.csv")
matrix260 <- read.csv(matrix260_path, check.names = FALSE)
stopifnot(colnames(matrix260)[1] == "ID", ncol(matrix260) == 13L)
rownames(matrix260) <- matrix260$ID
expr260 <- log2(as.matrix(matrix260[, -1, drop = FALSE]) + 1)
storage.mode(expr260) <- "numeric"
hape_cols <- grep("^hape", colnames(expr260))
control_cols <- grep("^con", colnames(expr260))
stopifnot(length(hape_cols) == 6L, length(control_cols) == 6L)
gene260 <- marker_panel
gene260$dataset <- "GSE260910"
gene260$comparison <- "descriptive HAPE-plus-batch versus control-plus-batch"
gene260$input_available <- gene260$gene %in% rownames(expr260)
gene260$expression_eligible <- NA
gene260$logFC <- vapply(gene260$gene, function(gene) {
  if (!gene %in% rownames(expr260)) return(NA_real_)
  mean(expr260[gene, hape_cols]) - mean(expr260[gene, control_cols])
}, numeric(1))
gene260$P.Value <- NA_real_
gene260$adj.P.Val <- NA_real_
gene260$direction <- direction_label(gene260$logFC,
                                     "higher_in_HAPE_plus_batch",
                                     "lower_in_HAPE_plus_batch")
gene260$significance <- "not estimable: disease and batch are completely confounded"
gene260$interpretation <- "descriptive group-plus-batch difference only; matrix unit unresolved"
gene260$whether_supports_composition_shift <- "no: disease-specific or composition-specific inference is prohibited"
gene260 <- gene260[, colnames(gene103)]

present260 <- lapply(marker_sets, function(x) intersect(x, rownames(expr260)))
scores260 <- make_scores(expr260, marker_sets)
score260 <- data.frame(
  dataset = "GSE260910",
  comparison = "descriptive HAPE-plus-batch versus control-plus-batch marker-score difference",
  cell_type = rownames(scores260$scores),
  marker_count = scores260$coverage$present_markers,
  effect = rowMeans(scores260$scores[, hape_cols, drop = FALSE]) -
    rowMeans(scores260$scores[, control_cols, drop = FALSE]),
  CI_lower = NA_real_, CI_upper = NA_real_, P_value = NA_real_, BH_FDR = NA_real_)
score260$direction <- direction_label(score260$effect,
                                      "higher_in_HAPE_plus_batch",
                                      "lower_in_HAPE_plus_batch")
score260$significance <- "not estimable: disease and batch are completely confounded"
score260$composition_interpretation <- "no HAPE-specific composition inference permitted"

gene_results <- rbind(gene103, gene103_b, gene756_avg, gene756_int, gene260)
score_results <- rbind(score103, score103_b, score756_avg, score756_int, score260)
write.csv(gene_results, file.path(out, "cell_type_marker_gene_results.csv"), row.names = FALSE)
write.csv(score_results, file.path(out, "cell_type_marker_score_results.csv"), row.names = FALSE)

coverage <- rbind(
  transform(scores103$coverage, dataset = "GSE103940"),
  transform(scores756$coverage, dataset = "GSE75665"),
  transform(scores260$coverage, dataset = "GSE260910")
)
coverage <- coverage[, c("dataset", setdiff(colnames(coverage), "dataset"))]
write.csv(coverage, file.path(out, "cell_type_marker_coverage.csv"), row.names = FALSE)

eligibility <- data.frame(
  dataset = c("GSE103940", "GSE75665", "GSE260910"),
  source = c("blood", "peripheral venous blood", "blood"),
  design = c("11 subjects; paired plain/high altitude",
             "10 subjects; repeated plain/high altitude; 5 AMS and 5 non-AMS",
             "12 independent subjects; 6 HAPE and 6 control"),
  composition_audit_eligibility = c(
    "eligible for paired marker-expression and marker-score audit",
    "eligible for repeated-measures marker-expression and marker-score audit",
    "descriptive only; disease and batch are completely confounded"),
  formal_deconvolution_performed = "No",
  reason_no_formal_deconvolution = c(
    "No matched CBC/flow-cytometry measurements or archived validated whole-blood signature/protocol; targeted marker audit avoids presenting proxy scores as proportions",
    "No matched CBC/flow-cytometry measurements or archived validated whole-blood signature/protocol; n=10 and the primary question is AMS-by-altitude interaction",
    "Complete disease-batch confounding prevents HAPE-specific interpretation regardless of deconvolution method")
)
write.csv(eligibility, file.path(out, "dataset_cell_composition_eligibility.csv"), row.names = FALSE)

# Supplementary candidate visualization: aggregate proxy effects only.
plot_df <- score_results[score_results$dataset != "GSE260910", ]
plot_marker_scores <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old))
  par(mfrow = c(4, 1), mar = c(4.2, 8.5, 2.5, 1), oma = c(4, 0, 0, 0),
      family = "sans")
  for (comparison in unique(plot_df$comparison)) {
    z <- plot_df[plot_df$comparison == comparison, ]
    z <- z[match(rev(names(marker_sets)), z$cell_type), ]
    limits <- range(c(z$CI_lower, z$CI_upper, 0), finite = TRUE)
    pad <- diff(limits) * 0.08
    if (!is.finite(pad) || pad == 0) pad <- 0.1
    plot(z$effect, seq_len(nrow(z)), xlim = limits + c(-pad, pad),
         ylim = c(0.5, nrow(z) + 0.5), yaxt = "n", ylab = "", xlab = "",
         pch = 19, cex = 1.1,
         col = ifelse(z$BH_FDR < 0.05, "#B33A3A", "#52606D"),
         main = comparison, cex.main = 0.9)
    axis(2, at = seq_len(nrow(z)), labels = z$cell_type, las = 1, cex.axis = 0.8)
    axis(1, cex.axis = 0.8)
    abline(v = 0, col = "grey75", lwd = 0.8)
    segments(z$CI_lower, seq_len(nrow(z)), z$CI_upper, seq_len(nrow(z)),
             col = ifelse(z$BH_FDR < 0.05, "#B33A3A", "#52606D"), lwd = 1.3)
  }
  mtext("Mean change in aggregate marker z-score (95% CI)", side = 1,
        outer = TRUE, line = 1.1, cex = 0.85)
  mtext("Expression-based lineage marker proxy; not an estimated cell proportion",
        side = 1, outer = TRUE, line = 2.5, cex = 0.72)
}
pdf(file.path(out, "Supplementary_cell_composition_marker_score_audit.pdf"),
    width = 7.1, height = 10.8, family = "Helvetica")
plot_marker_scores()
dev.off()
png(file.path(out, "Supplementary_cell_composition_marker_score_audit.png"),
    width = 4260, height = 6480, res = 600, type = "cairo")
plot_marker_scores()
dev.off()

writeLines(capture.output(sessionInfo()), file.path(out, "cell_composition_audit_sessionInfo.txt"))

cat("GSE103940 marker-score results\n")
print(score103)
cat("\nGSE75665 average altitude marker-score results\n")
print(score756_avg)
cat("\nGSE75665 interaction marker-score results\n")
print(score756_int)
