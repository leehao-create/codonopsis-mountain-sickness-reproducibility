# GSE260910 disease-batch confounding report

## Decision

**Conservative classification: Case C, complete disease-batch confounding. Formal DEG analysis was stopped.**

GEO provides disease state, GSM, original sample description and platform for all 12 samples, but no sample-level sequencing-batch field. The source publication states that the HAPE and control groups corresponded to two sequencing batches. Exact physical batch identifiers were not deposited. Because the only published batch relation is group-determining, disease and batch cannot be entered as independently estimable effects. No revised HAPE DEG table was generated.

Source records:

- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE260910
- Publication: https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1444666/full
- Frozen local GEO SOFT: `source_metadata/GSE260910_family.soft.gz`

## Sample/group/batch distribution

| Reported sequencing-batch relation | HAPE | control |
|---|---:|---:|
| batch corresponding to HAPE | 6 | 0 |
| batch corresponding to control | 0 | 6 |

The exact batch names are deliberately not invented. The labels above encode only the published group-to-batch correspondence.

## Design-matrix rank

For `~ disease_status + sequencing_batch`, the model matrix has 3 columns but rank 2.
The disease indicator and batch indicator are identical, so one is an exact linear combination of the other. A disease coefficient conditional on batch does not exist as an independently estimable parameter.

## Why batch correction is prohibited

ComBat, `removeBatchEffect()`, Harmony, or adding batch to the design cannot determine which portion of the same contrast is biological and which is technical. With no within-batch HAPE/control comparison, correction would either remove the disease contrast or assign an arbitrary interpretation to it.

## Affected manuscript evidence

- Existing HAPE DEG counts and directions cannot be retained as disease-specific confirmatory evidence.
- HAPE GO/KEGG and other pathway results derived from those DEGs are affected.
- The HAPE PPI/Cytoscape network is affected because its node list and directions derive from those DEGs.
- HAPE candidate-expression comparisons, including MMP9, can reflect batch.
- GSE260910 must not enter supervised machine-learning feature selection or performance evaluation.
- The bulk-derived panels of original Figure 4 must be removed, replaced with a valid independent HAPE dataset, or explicitly reframed as non-inferential/confounded exploratory material. The lung reference-atlas panels are a separate provenance question.

## Evidence boundary

The publication-level wording is sufficient to prohibit an unqualified disease-effect analysis, but the public records do not expose physical batch IDs per sample. If a verified laboratory batch sheet later demonstrates that both HAPE and controls occurred within every batch, this decision can be revisited from the frozen metadata; until then, no HAPE DEG is statistically defensible.
