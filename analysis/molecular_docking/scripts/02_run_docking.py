#!/usr/bin/env python3
"""Prepare ligands, validate pockets by redocking, and dock four fixed pairs."""

from pathlib import Path
import csv
import subprocess
import shutil

from rdkit import Chem
from rdkit.Chem import rdMolAlign


ROOT = Path(__file__).resolve().parents[1]
PREP = ROOT / "prepared"
RESULTS = ROOT / "results"
CONFIGS = ROOT / "configs"
LOGS = ROOT / "logs"
for directory in (RESULTS, CONFIGS, LOGS):
    directory.mkdir(parents=True, exist_ok=True)

SEED = 20260927
NUM_MODES = 20
ENERGY_RANGE = 8
REDOCK_EXHAUSTIVENESS = 64
CANDIDATE_EXHAUSTIVENESS = 32

PAIRS = [
    {"target": "NR3C2", "ligand_id": "MOL000449", "ligand": "stigmasterol", "pubchem_cid": 5280794},
    {"target": "NR3C2", "ligand_id": "MOL004355", "ligand": "spinasterol", "pubchem_cid": 5281331},
    {"target": "NR3C2", "ligand_id": "MOL008407", "ligand": "steroid-like TCMSP compound", "pubchem_cid": 14807783},
    {"target": "MMP9", "ligand_id": "MOL000006", "ligand": "luteolin", "pubchem_cid": 5280445},
]


def run(cmd, log_path):
    result = subprocess.run(cmd, text=True, capture_output=True)
    log_path.write_text("COMMAND: " + " ".join(map(str, cmd)) + "\n\n" + result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f"Command failed: {' '.join(map(str, cmd))}\n{result.stderr}")


def write_config(path, receptor, ligand, output, box, exhaustiveness):
    values = {
        "receptor": receptor, "ligand": ligand, "out": output,
        "center_x": box["center_x"], "center_y": box["center_y"], "center_z": box["center_z"],
        "size_x": box["size_x"], "size_y": box["size_y"], "size_z": box["size_z"],
        "exhaustiveness": exhaustiveness, "num_modes": NUM_MODES,
        "energy_range": ENERGY_RANGE, "seed": SEED, "cpu": 8,
    }
    path.write_text("\n".join(f"{key} = {value}" for key, value in values.items()) + "\n")


def scores(path):
    return [float(line.split()[3]) for line in path.read_text().splitlines()
            if line.startswith("REMARK VINA RESULT:")]


with (PREP / "receptor_and_grid_metadata.csv").open() as handle:
    boxes = {row["target"]: row for row in csv.DictReader(handle)}

redocking_rows = []
for target, box in boxes.items():
    outdir = RESULTS / "redocking" / target
    outdir.mkdir(parents=True, exist_ok=True)
    native_protonated = outdir / "native_ligand_pH7.4.sdf"
    native_pdbqt = outdir / "native_ligand.pdbqt"
    native_out = outdir / "native_redocked.pdbqt"
    native_sdf = outdir / "native_redocked.sdf"

    run([
        shutil.which("obabel") or "obabel", "-isdf", str(PREP / target / "native_ligand_crystal.sdf"),
        "-osdf", "-O", str(native_protonated), "-p", "7.4",
    ], LOGS / f"02_{target}_native_protonation.log")
    run([
        shutil.which("mk_prepare_ligand.py") or "mk_prepare_ligand.py", "-i", str(native_protonated),
        "-o", str(native_pdbqt), "--charge_model", "gasteiger", "--add_index_map",
    ], LOGS / f"02_{target}_native_meeko.log")
    config = CONFIGS / f"redock_{target}.txt"
    write_config(config, PREP / target / "receptor.pdbqt", native_pdbqt, native_out,
                 box, REDOCK_EXHAUSTIVENESS)
    run([shutil.which("vina") or "vina", "--config", str(config)], LOGS / f"02_{target}_redocking_vina.log")
    run([shutil.which("mk_export.py") or "mk_export.py", str(native_out), "-s", str(native_sdf)],
        LOGS / f"02_{target}_redocking_export.log")

    reference = Chem.RemoveHs(Chem.SDMolSupplier(
        str(PREP / target / "native_ligand_crystal.sdf"), removeHs=False)[0])
    poses = [Chem.RemoveHs(mol) for mol in Chem.SDMolSupplier(str(native_sdf), removeHs=False)
             if mol is not None]
    affinities = scores(native_out)
    rmsds = [float(rdMolAlign.GetBestRMS(reference, pose)) for pose in poses]
    if len(rmsds) != len(affinities) or not rmsds:
        raise RuntimeError(f"Redocking parse mismatch for {target}")
    with (outdir / "all_pose_rmsd_affinity.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["pose_rank", "affinity_kcal_mol", "RMSD_A"])
        writer.writeheader()
        writer.writerows({"pose_rank": i + 1, "affinity_kcal_mol": affinity, "RMSD_A": rmsd}
                         for i, (affinity, rmsd) in enumerate(zip(affinities, rmsds)))
    redocking_rows.append({
        "target": target, "pdb_id": box["pdb_id"],
        "co_crystallized_ligand": box["co_crystallized_ligand"],
        "top_pose_affinity_kcal_mol": affinities[0],
        "top_pose_symmetry_corrected_heavy_atom_RMSD_A": round(rmsds[0], 3),
        "minimum_RMSD_A": round(min(rmsds), 3), "minimum_RMSD_pose_rank": rmsds.index(min(rmsds)) + 1,
        "pass_threshold_A": 2.0, "top_pose_pass": rmsds[0] <= 2.0,
        "RMSD_method": "RDKit GetBestRMS, symmetry-corrected heavy-atom RMSD in unchanged receptor frame",
        "exhaustiveness": REDOCK_EXHAUSTIVENESS, "num_modes": NUM_MODES,
        "energy_range": ENERGY_RANGE, "random_seed": SEED,
    })

with (ROOT / "redocking_validation.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=redocking_rows[0].keys())
    writer.writeheader()
    writer.writerows(redocking_rows)

candidate_rows = []
for pair in PAIRS:
    target, ligand_id = pair["target"], pair["ligand_id"]
    box = boxes[target]
    outdir = RESULTS / "candidate_docking" / target / ligand_id
    outdir.mkdir(parents=True, exist_ok=True)
    source = ROOT / "inputs" / "candidate_ligands" / f"{ligand_id}_PubChem3D.sdf"
    minimized = outdir / f"{ligand_id}_pH7.4_MMFF94_minimized.sdf"
    ligand_pdbqt = outdir / f"{ligand_id}.pdbqt"
    docked_pdbqt = outdir / f"{ligand_id}_docked.pdbqt"
    docked_sdf = outdir / f"{ligand_id}_docked.sdf"

    run([
        shutil.which("obabel") or "obabel", "-isdf", str(source), "-osdf", "-O", str(minimized),
        "-p", "7.4", "--minimize", "--ff", "MMFF94", "--steps", "500", "--crit", "1e-6",
    ], LOGS / f"03_{target}_{ligand_id}_protonation_minimization.log")
    run([
        shutil.which("mk_prepare_ligand.py") or "mk_prepare_ligand.py", "-i", str(minimized), "-o", str(ligand_pdbqt),
        "--charge_model", "gasteiger", "--add_index_map",
    ], LOGS / f"03_{target}_{ligand_id}_meeko.log")
    config = CONFIGS / f"dock_{target}_{ligand_id}.txt"
    write_config(config, PREP / target / "receptor.pdbqt", ligand_pdbqt, docked_pdbqt,
                 box, CANDIDATE_EXHAUSTIVENESS)
    run([shutil.which("vina") or "vina", "--config", str(config)],
        LOGS / f"03_{target}_{ligand_id}_vina.log")
    run([shutil.which("mk_export.py") or "mk_export.py", str(docked_pdbqt), "-s", str(docked_sdf)],
        LOGS / f"03_{target}_{ligand_id}_export.log")
    affinities = scores(docked_pdbqt)
    if not affinities:
        raise RuntimeError(f"No affinity parsed for {target}-{ligand_id}")
    candidate_rows.append({
        **pair, "ligand_source": f"archived PubChem 3D SDF, CID {pair['pubchem_cid']}",
        "ligand_preparation": "Open Babel pH 7.4; MMFF94, 500 steps, convergence 1e-6; Meeko Gasteiger charges",
        "pdb_id": box["pdb_id"], "chain": box["chain"],
        "binding_site": f"co-crystallized {box['co_crystallized_ligand']} pocket",
        "grid_center_xyz": f"{box['center_x']};{box['center_y']};{box['center_z']}",
        "grid_size_xyz_A": f"{box['size_x']};{box['size_y']};{box['size_z']}",
        "vina_version": "conda package 1.2.7; executable reports AutoDock Vina f458505-mod",
        "exhaustiveness": CANDIDATE_EXHAUSTIVENESS, "num_modes": NUM_MODES,
        "energy_range": ENERGY_RANGE, "random_seed": SEED,
        "modes_generated": len(affinities), "best_affinity_kcal_mol": affinities[0],
        "pose_selection_rule": "lowest Vina score (rank 1); contacts reported descriptively",
        "best_pose_file": str(docked_sdf),
    })

with (ROOT / "new_docking_results.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=candidate_rows[0].keys())
    writer.writeheader()
    writer.writerows(candidate_rows)
