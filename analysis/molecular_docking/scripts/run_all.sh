#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$root/logs"
cd "$root"

bash "$root/scripts/00_fetch_inputs.sh" > "$root/logs/00_fetch_inputs.log" 2>&1
python "$root/scripts/01_prepare_receptors.py" > "$root/logs/01_prepare_receptors.log" 2>&1
python "$root/scripts/02_run_docking.py" > "$root/logs/02_run_docking.log" 2>&1
python "$root/scripts/03_analyze_and_compare.py" > "$root/logs/03_analyze_and_compare.log" 2>&1
python "$root/scripts/04_validate_outputs.py" > "$root/logs/04_validate_outputs.log" 2>&1

{
  python --version
  vina --version
  obabel -V
  python -c 'import importlib.metadata as m, meeko, rdkit, openmm; print("Meeko", meeko.__version__); print("RDKit", rdkit.__version__); print("OpenMM", openmm.__version__); print("PDBFixer", m.version("pdbfixer"))'
} > "$root/software_versions.txt" 2>&1

(cd "$root" && find . -type f ! -name output_manifest_sha256.txt -print0 | sort -z | xargs -0 sha256sum) \
  > "$root/output_manifest_sha256.txt"
