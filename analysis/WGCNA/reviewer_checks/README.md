# Reviewer-requested WGCNA checks

These analyses were performed as reviewer-requested sensitivity or diagnostic checks and were not used to redefine the primary manuscript WGCNA analysis.

This directory contains the binary-AMS 12-test check, the original 12-module by four-group-indicator (48-test) matrices and BH adjustment, module-size-aware candidate-target enrichment, and post hoc MMP9 diagnostics. Alternative network-type and soft-threshold investigations belong here rather than under `manuscript_primary/`.

`scripts/02_MMP9_diagnostics_from_frozen_modules.R` documents the calculation of traditional GS/P, kME and signed-adjacency kWithin from the frozen expression matrix and frozen module assignments. It does not reconstruct or alter modules.
