options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
module_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
out <- Sys.getenv("DE_OUTPUT_DIR", file.path(module_dir, "reproduced_results", "GSE260910_confounding"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
map_file <- file.path(module_dir, "data", "metadata", "GSE260910_sample_map.csv")
soft_file <- Sys.getenv("GSE260910_SOFT", file.path(module_dir, "data", "external", "GSE260910_family.soft.gz"))
stopifnot(dir.exists(out), file.exists(map_file), file.exists(soft_file))

sample_map <- read.csv(map_file, check.names = FALSE)
soft <- readLines(gzfile(soft_file), warn = FALSE)
stopifnot(nrow(sample_map) == 12L, !anyDuplicated(sample_map$gsm),
          all(sample_map$gsm %in% sub("^\\^SAMPLE = ", "", grep("^\\^SAMPLE = ", soft, value = TRUE))),
          !any(grepl("batch", soft, ignore.case = TRUE)))

# GEO supplies no sample-level batch field. The source paper states that HAPE and
# control groups corresponded to two sequencing batches. We encode only that
# published correspondence, without inventing physical batch identifiers.
sample_map$sequencing_batch <- ifelse(
  sample_map$disease_status == "HAPE",
  "reported_batch_corresponding_to_HAPE",
  "reported_batch_corresponding_to_control"
)
sample_map$batch_source <- paste(
  "Publication-level group-to-batch correspondence;",
  "exact sample-level physical batch IDs are not deposited in GEO"
)
sample_map$batch_assignment_confidence <- "group-level statement; exact batch ID unavailable"
write.csv(sample_map, file.path(out, "GSE260910_sample_disease_batch_table.csv"), row.names = FALSE)

tab <- as.data.frame.matrix(table(sample_map$sequencing_batch, sample_map$disease_status))
tab$sequencing_batch <- rownames(tab)
tab <- tab[, c("sequencing_batch", setdiff(colnames(tab), "sequencing_batch"))]
write.csv(tab, file.path(out, "GSE260910_disease_batch_contingency.csv"), row.names = FALSE)

model_data <- data.frame(
  disease_status = relevel(factor(sample_map$disease_status), ref = "control"),
  sequencing_batch = relevel(factor(sample_map$sequencing_batch),
                              ref = "reported_batch_corresponding_to_control")
)
design <- model.matrix(~ disease_status + sequencing_batch, model_data)
rank <- qr(design)$rank
rank_summary <- data.frame(
  samples = nrow(design), design_columns = ncol(design), design_rank = rank,
  full_column_rank = rank == ncol(design),
  disease_batch_complete_confounding = TRUE,
  decision = "Case C: stop formal GSE260910 DEG"
)
write.csv(data.frame(sample_id = sample_map$gsm, design, check.names = FALSE),
          file.path(out, "GSE260910_proposed_design_matrix_rank_deficient.csv"), row.names = FALSE)
write.csv(rank_summary, file.path(out, "GSE260910_design_rank_summary.csv"), row.names = FALSE)

affected <- data.frame(
  analysis = c("HAPE DEG", "HAPE pathway enrichment", "HAPE PPI/network",
               "HAPE candidate-gene expression", "HAPE supervised ML", "Manuscript Figure 4"),
  existing_source = c(
    "lhh/04.bulk/2/hape.con.DEG.csv",
    "lhh/04.bulk/2/hape.con.DEG.david.xlsx; david.*.txt; dot.pdf",
    "lhh/04.bulk/2/network.csv; up.down.cytoscape.tsv; cyc.cys",
    "lhh/04.bulk/2/logcpm.matrix.xlsx",
    "lhh/07.机器学习 outputs/history",
    "lhh/11.投稿/figure4.pdf (HAPE panels)"
  ),
  consequence = c(
    "Disease effect cannot be separated from reported batch effect; do not use as confirmatory HAPE DEG evidence",
    "Derived from non-identifiable DEG signal; enrichment is not disease-specific evidence",
    "Derived from affected DEGs; network interpretation is not disease-specific",
    "Group differences may be batch-driven",
    "Prohibited: classifier may learn sequencing batch rather than HAPE",
    "Bulk HAPE DEG, pathway, PPI and candidate-expression panels require removal or explicit non-inferential limitation"
  )
)
write.csv(affected, file.path(out, "GSE260910_affected_outputs.csv"), row.names = FALSE)

report <- c(
  "# GSE260910 disease-batch confounding report",
  "",
  "## Decision",
  "",
  "**Conservative classification: Case C, complete disease-batch confounding. Formal DEG analysis was stopped.**",
  "",
  "GEO provides disease state, GSM, original sample description and platform for all 12 samples, but no sample-level sequencing-batch field. The source publication states that the HAPE and control groups corresponded to two sequencing batches. Exact physical batch identifiers were not deposited. Because the only published batch relation is group-determining, disease and batch cannot be entered as independently estimable effects. No revised HAPE DEG table was generated.",
  "",
  "Source records:",
  "",
  "- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE260910",
  "- Publication: https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1444666/full",
  "- GEO family SOFT supplied through `GSE260910_SOFT` (not redistributed)",
  "",
  "## Sample/group/batch distribution",
  "",
  "| Reported sequencing-batch relation | HAPE | control |",
  "|---|---:|---:|",
  "| batch corresponding to HAPE | 6 | 0 |",
  "| batch corresponding to control | 0 | 6 |",
  "",
  "The exact batch names are deliberately not invented. The labels above encode only the published group-to-batch correspondence.",
  "",
  "## Design-matrix rank",
  "",
  paste0("For `~ disease_status + sequencing_batch`, the model matrix has 3 columns but rank ", rank, "."),
  "The disease indicator and batch indicator are identical, so one is an exact linear combination of the other. A disease coefficient conditional on batch does not exist as an independently estimable parameter.",
  "",
  "## Why batch correction is prohibited",
  "",
  "ComBat, `removeBatchEffect()`, Harmony, or adding batch to the design cannot determine which portion of the same contrast is biological and which is technical. With no within-batch HAPE/control comparison, correction would either remove the disease contrast or assign an arbitrary interpretation to it.",
  "",
  "## Affected manuscript evidence",
  "",
  "- Existing HAPE DEG counts and directions cannot be retained as disease-specific confirmatory evidence.",
  "- HAPE GO/KEGG and other pathway results derived from those DEGs are affected.",
  "- The HAPE PPI/Cytoscape network is affected because its node list and directions derive from those DEGs.",
  "- HAPE candidate-expression comparisons, including MMP9, can reflect batch.",
  "- GSE260910 must not enter supervised machine-learning feature selection or performance evaluation.",
  "- The bulk-derived panels of original Figure 4 must be removed, replaced with a valid independent HAPE dataset, or explicitly reframed as non-inferential/confounded exploratory material. The lung reference-atlas panels are a separate provenance question.",
  "",
  "## Evidence boundary",
  "",
  "The publication-level wording is sufficient to prohibit an unqualified disease-effect analysis, but the public records do not expose physical batch IDs per sample. If a verified laboratory batch sheet later demonstrates that both HAPE and controls occurred within every batch, this decision can be revisited from the frozen metadata; until then, no HAPE DEG is statistically defensible."
)
writeLines(report, file.path(out, "GSE260910_batch_confounding_report.md"), useBytes = TRUE)

write.csv(data.frame(
  path = c(map_file, soft_file),
  md5 = unname(tools::md5sum(c(map_file, soft_file)))
), file.path(out, "GSE260910_audit_input_md5.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "GSE260910_batch_audit_sessionInfo.txt"))
message("GSE260910 classified as Case C; formal DEG stopped. Design rank ", rank,
        "/", ncol(design), ".")
