# Configuration-file note

The CA2/BCL2 workflow did not generate standalone Vina `.txt` configuration files. It passed the fixed parameters directly from `scripts/02_run_redocking.py` and `scripts/03_run_candidate_docking.py`. The receptor-specific grid centers and sizes are stored in `../prepared/docking_boxes.csv`, while the seed, exhaustiveness, number of modes and energy range are preserved in the scripts and frozen Vina logs.
