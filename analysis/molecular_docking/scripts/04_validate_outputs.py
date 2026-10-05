#!/usr/bin/env python3
"""Fail closed if any required docking output is incomplete or invalid."""

from pathlib import Path
import csv


ROOT = Path(__file__).resolve().parents[1]

with (ROOT / "redocking_validation.csv").open() as handle:
    validation = list(csv.DictReader(handle))
assert len(validation) == 2
assert all(row["top_pose_pass"] == "True" for row in validation)
assert all(float(row["top_pose_symmetry_corrected_heavy_atom_RMSD_A"]) <= 2.0 for row in validation)

with (ROOT / "new_docking_results.csv").open() as handle:
    results = list(csv.DictReader(handle))
assert len(results) == 4
assert {(row["target"], row["ligand_id"]) for row in results} == {
    ("NR3C2", "MOL000449"), ("NR3C2", "MOL004355"),
    ("NR3C2", "MOL008407"), ("MMP9", "MOL000006"),
}
assert all(float(row["best_affinity_kcal_mol"]) != 0.0 for row in results)
assert all((ROOT / row["best_pose_file"]).is_file() for row in results)

for target in ("NR3C2", "MMP9"):
    receptor = ROOT / "prepared" / target / "receptor.pdbqt"
    assert receptor.stat().st_size > 10000
for config in ROOT.glob("configs/*.txt"):
    text = config.read_text()
    for field in ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z",
                  "exhaustiveness", "num_modes", "energy_range", "seed"):
        assert f"{field} =" in text

print("PASS: 2 redocking validations, 4 candidate pairs, nonempty receptors, complete configs")
