#!/usr/bin/env python3
"""Prepare native ligands, run fixed-seed Vina redocking and compute RMSD."""

from pathlib import Path
import csv
import re
import subprocess

from rdkit import Chem
from rdkit.Chem import rdMolAlign


ROOT = Path(__file__).resolve().parents[1]
PREP = ROOT / "prepared"
RESULTS = ROOT / "results" / "redocking"
LOGS = ROOT / "logs"
RESULTS.mkdir(parents=True, exist_ok=True)
LOGS.mkdir(exist_ok=True)

SEED = 20260917
EXHAUSTIVENESS = 64
NUM_MODES = 20
RMSD_PASS = 2.0


def run(cmd, log_file=None):
    result = subprocess.run(cmd, text=True, capture_output=True)
    if log_file:
        log_file.write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(map(str, cmd))}\n{result.stderr}")
    return result


def affinities_from_pdbqt(path):
    vals = []
    for line in path.read_text().splitlines():
        if line.startswith("REMARK VINA RESULT:"):
            vals.append(float(line.split()[3]))
    return vals


with (PREP / "docking_boxes.csv").open() as handle:
    boxes = {row["target"]: row for row in csv.DictReader(handle)}

rows = []
for target, box in boxes.items():
    tprep = PREP / target
    outdir = RESULTS / target
    outdir.mkdir(exist_ok=True)

    protonated = tprep / "native_ligand_pH7.4.sdf"
    run([
        "obabel", "-isdf", str(tprep / "native_ligand_crystal.sdf"),
        "-osdf", "-O", str(protonated), "-p", "7.4",
    ], LOGS / f"02_{target}_native_protonation.log")

    ligand_pdbqt = tprep / "native_ligand.pdbqt"
    run([
        "mk_prepare_ligand.py", "-i", str(protonated),
        "-o", str(ligand_pdbqt), "--charge_model", "gasteiger", "--add_index_map",
    ], LOGS / f"02_{target}_native_ligand_prep.log")

    docked_pdbqt = outdir / "native_redocked.pdbqt"
    vina_cmd = [
        "vina", "--receptor", str(tprep / "receptor.pdbqt"),
        "--ligand", str(ligand_pdbqt), "--out", str(docked_pdbqt),
        "--center_x", box["center_x"], "--center_y", box["center_y"],
        "--center_z", box["center_z"], "--size_x", box["size_x"],
        "--size_y", box["size_y"], "--size_z", box["size_z"],
        "--exhaustiveness", str(EXHAUSTIVENESS), "--seed", str(SEED),
        "--num_modes", str(NUM_MODES), "--energy_range", "8", "--cpu", "8",
    ]
    run(vina_cmd, LOGS / f"02_{target}_redocking_vina.log")

    docked_sdf = outdir / "native_redocked.sdf"
    run([
        "mk_export.py", str(docked_pdbqt), "-s", str(docked_sdf),
    ], LOGS / f"02_{target}_redocking_export.log")

    reference = Chem.RemoveHs(Chem.SDMolSupplier(str(tprep / "native_ligand_crystal.sdf"), removeHs=False)[0])
    docked = [Chem.RemoveHs(m) for m in Chem.SDMolSupplier(str(docked_sdf), removeHs=False) if m is not None]
    rmsds = [float(rdMolAlign.GetBestRMS(reference, pose)) for pose in docked]
    affinities = affinities_from_pdbqt(docked_pdbqt)
    if not rmsds or not affinities:
        raise RuntimeError(f"No redocking poses parsed for {target}")
    if len(affinities) != len(rmsds):
        raise RuntimeError(f"Pose/affinity count mismatch for {target}: {len(rmsds)} vs {len(affinities)}")

    with (outdir / "pose_rmsd_and_affinity.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["pose_rank", "affinity_kcal_mol", "symmetry_corrected_heavy_atom_RMSD_A"])
        writer.writeheader()
        writer.writerows({"pose_rank": i + 1, "affinity_kcal_mol": a, "symmetry_corrected_heavy_atom_RMSD_A": r}
                         for i, (a, r) in enumerate(zip(affinities, rmsds)))

    rows.append({
        "target": target,
        "pdb_id": box["pdb_id"],
        "chain": box["chain"],
        "co_crystallized_ligand": box["native_ligand"],
        "observed_top_pose_affinity_kcal_mol": affinities[0],
        "top_ranked_pose_RMSD_A": round(rmsds[0], 3),
        "minimum_RMSD_among_generated_poses_A": round(min(rmsds), 3),
        "minimum_RMSD_pose_rank": rmsds.index(min(rmsds)) + 1,
        "validation_threshold_A": RMSD_PASS,
        "protocol_pass": rmsds[0] <= RMSD_PASS,
        "rmsd_method": "RDKit symmetry-corrected heavy-atom GetBestRMS; receptor frame unchanged",
        "protonation_pH": 7.4,
        "charge_model": "Gasteiger (Meeko 0.8.0)",
        "exhaustiveness": EXHAUSTIVENESS,
        "seed": SEED,
        "num_modes": NUM_MODES,
        "grid_center": f"{box['center_x']};{box['center_y']};{box['center_z']}",
        "grid_size_A": f"{box['size_x']};{box['size_y']};{box['size_z']}",
    })

with (ROOT / "redocking_validation.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
