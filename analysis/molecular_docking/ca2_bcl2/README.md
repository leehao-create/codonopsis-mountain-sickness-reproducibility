# CA2 and BCL2 docking archive

This directory contains the frozen files actually used for the CA2 (PDB 5E28, chain A) and BCL2 (PDB 8HTS, chain A) results reported in the manuscript. No docking was rerun while assembling this release.

## Contents

- `inputs/`: RCSB receptor/native-ligand records and the archived PubChem/TCMSP ligand structures.
- `prepared/`: receptor and ligand preparation outputs. `prepared/docking_boxes.csv` is the authoritative pocket/grid configuration record.
- `results/redocking/`: all redocked poses and per-pose RMSD/affinity tables.
- `results/candidate_docking/`: all generated candidate poses in PDBQT and SDF formats.
- `logs/`: frozen receptor/ligand preparation and Vina execution logs.
- `scripts/`: the original analysis logic with only local absolute paths replaced by repository-relative paths and active-environment commands.
- `redocking_validation.csv`: frozen top-pose validation metrics and parameters.
- `new_docking_results.csv` and `new_docking_interactions.csv`: frozen candidate scores, pose paths and geometric contact summaries.
- `source_provenance.csv`: receptor, ligand and compound-target mapping provenance.
- `figures/manuscript_Fig6_final.pdf`: the archived manuscript Figure 6 used to trace the reported panels to these results.

## Fixed configuration

CA2 used the BC5-defined pocket centered at `(-3.770, 4.124, 14.906)` A with a `22 x 22 x 22` A box. BCL2 used the N2L-defined pocket centered at `(-0.210, 9.374, 7.186)` A with a `26 x 27 x 27` A box. AutoDock Vina 1.2.7 used seed `20260917`, `num_modes=20`, `energy_range=8`, exhaustiveness 64 for redocking and 32 for candidate docking.

The original workflow invoked Vina through Python/subprocess commands and did not create standalone Vina configuration text files. The exact equivalent settings are retained in `prepared/docking_boxes.csv`, `scripts/02_run_redocking.py`, `scripts/03_run_candidate_docking.py`, and the Vina logs. These files are the configuration record; no retrospective config file was fabricated.

The archived CA2/BCL2 scripts and logs document PubChem 3D input, pH 7.4 protonation and Gasteiger charge/PDBQT preparation. They do not contain a separate MMFF94 energy-minimization command or log, so no such missing record was inferred.

## Frozen outcomes

- CA2/BC5 top-pose redocking RMSD: 0.504 A.
- BCL2/N2L top-pose redocking RMSD: 0.339 A.
- CA2/MOL007059 top predicted affinity: -8.706 kcal/mol.
- CA2/MOL008413 top predicted affinity: -5.214 kcal/mol.
- BCL2/MOL002579 (capsaicin) top predicted affinity: -7.271 kcal/mol.
- BCL2/MOL000008 (apigenin) top predicted affinity: -6.924 kcal/mol.

These are frozen computational pose-compatibility outputs, not evidence of physical binding or biological activity.
