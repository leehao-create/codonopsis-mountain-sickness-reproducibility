# Molecular docking

This module contains the frozen reproducibility materials for all four receptor structures reported in the manuscript:

- NR3C2: PDB 2AA2, chain A.
- MMP9: PDB 1GKC, chain A.
- CA2: PDB 5E28, chain A.
- BCL2: PDB 8HTS, chain A.

The module root contains the NR3C2/MMP9 receptor and ligand inputs, prepared structures, standalone Vina configurations, execution logs, candidate poses, redocking poses, RMSD results and scripts. Create the environment from `environment.yml`, activate it, and run `scripts/run_all.sh` only if reproduction is required.

The complete CA2/BCL2 archive is isolated under `ca2_bcl2/`. It includes inputs, prepared receptors/ligands, all candidate and redocking poses, per-pose RMSD tables, logs, fixed parameters, provenance and the archived manuscript Figure 6. Its workflow passed parameters directly through the Python commands, so no standalone Vina `.txt` configuration existed; `ca2_bcl2/prepared/docking_boxes.csv`, the two docking scripts and Vina logs jointly preserve the exact configuration.

CA2/BCL2 Vina parameters are recorded in the analysis scripts, `docking_boxes.csv`, and corresponding docking logs. Ligand preparation procedures are documented in the scripts and provenance records.

`compound_target_provenance.csv` links all reported receptor-compound pairs to their archived project mapping source, PubChem identity and public frozen result record. Raw TCMSP mapping exports are not redistributed.

No docking calculation was rerun during repository assembly. Numerical results remain frozen in the CSV and pose files.
