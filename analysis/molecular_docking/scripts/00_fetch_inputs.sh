#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$root/inputs/structures" "$root/inputs/native_ligands" \
  "$root/inputs/candidate_ligands" "$root/logs"

curl --fail --location --retry 3 --silent --show-error \
  "https://files.rcsb.org/download/2AA2.pdb" \
  --output "$root/inputs/structures/2AA2.pdb"
curl --fail --location --retry 3 --silent --show-error \
  "https://files.rcsb.org/download/1GKC.pdb" \
  --output "$root/inputs/structures/1GKC.pdb"
curl --fail --location --retry 3 --silent --show-error \
  "https://files.rcsb.org/ligands/download/AS4_ideal.sdf" \
  --output "$root/inputs/native_ligands/AS4_ideal.sdf"
curl --fail --location --retry 3 --silent --show-error \
  "https://files.rcsb.org/ligands/download/NFH_ideal.sdf" \
  --output "$root/inputs/native_ligands/NFH_ideal.sdf"

for mol in MOL000449 MOL004355 MOL008407 MOL000006; do
  test -s "$root/inputs/candidate_ligands/${mol}_PubChem3D.sdf"
done

sha256sum "$root"/inputs/structures/* "$root"/inputs/native_ligands/* \
  "$root"/inputs/candidate_ligands/* > "$root/inputs/input_sha256.txt"
