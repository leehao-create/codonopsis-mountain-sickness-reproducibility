#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$root/inputs/structures" "$root/inputs/ligands" "$root/logs"

curl --fail --location --retry 3 --silent --show-error \
  "https://files.rcsb.org/download/5E28.pdb" \
  --output "$root/inputs/structures/5E28.pdb"
curl --fail --location --retry 3 --silent --show-error \
  "https://files.rcsb.org/download/8HTS.pdb" \
  --output "$root/inputs/structures/8HTS.pdb"

for component in BC5 N2L; do
  curl --fail --location --retry 3 --silent --show-error \
    "https://files.rcsb.org/ligands/download/${component}_ideal.sdf" \
    --output "$root/inputs/ligands/${component}_ideal.sdf"
done

for mol in MOL007059 MOL008413 MOL002579 MOL000008; do
  curl --fail --location --retry 3 --silent --show-error \
    "https://www.tcmsp-e.com/tcmspmol/${mol}.mol2" \
    --output "$root/inputs/ligands/${mol}_TCMSP.mol2"
done

declare -A pubchem=(
  [MOL007059]=5318290
  [MOL008413]=178371
  [MOL002579]=1548943
  [MOL000008]=5280443
)
for mol in "${!pubchem[@]}"; do
  cid="${pubchem[$mol]}"
  curl --fail --location --retry 3 --silent --show-error \
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/SDF?record_type=3d" \
    --output "$root/inputs/ligands/${mol}_PubChem3D.sdf"
done

sha256sum "$root"/inputs/structures/* "$root"/inputs/ligands/* \
  > "$root/inputs/input_sha256.txt"
