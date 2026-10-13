#!/usr/bin/env python3
"""Dock only mapped CA2/BCL2 candidates after successful redocking."""

from pathlib import Path
import csv
import subprocess


ROOT = Path(__file__).resolve().parents[1]
PREP = ROOT / "prepared"
RESULTS = ROOT / "results" / "candidate_docking"
LOGS = ROOT / "logs"
RESULTS.mkdir(parents=True, exist_ok=True)

SEED = 20260917
EXHAUSTIVENESS = 32
NUM_MODES = 20

CANDIDATES = [
    {"target": "CA2", "ligand_id": "MOL007059", "ligand": "3-beta-Hydroxymethyllenetanshiquinone", "pubchem_cid": 5318290},
    {"target": "CA2", "ligand_id": "MOL008413", "ligand": "Codopiloic acid", "pubchem_cid": 178371},
    {"target": "BCL2", "ligand_id": "MOL002579", "ligand": "capsaicin", "pubchem_cid": 1548943},
    {"target": "BCL2", "ligand_id": "MOL000008", "ligand": "apigenin", "pubchem_cid": 5280443},
]


def run(cmd, log_file):
    result = subprocess.run(cmd, text=True, capture_output=True)
    log_file.write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(map(str, cmd))}\n{result.stderr}")


def affinities(path):
    return [float(line.split()[3]) for line in path.read_text().splitlines()
            if line.startswith("REMARK VINA RESULT:")]


with (ROOT / "redocking_validation.csv").open() as handle:
    validation = {row["target"]: row["protocol_pass"].lower() == "true" for row in csv.DictReader(handle)}
with (PREP / "docking_boxes.csv").open() as handle:
    boxes = {row["target"]: row for row in csv.DictReader(handle)}

rows = []
for item in CANDIDATES:
    target = item["target"]
    if not validation.get(target, False):
        raise RuntimeError(f"Stopping rule triggered: redocking did not pass for {target}")
    box = boxes[target]
    outdir = RESULTS / target / item["ligand_id"]
    outdir.mkdir(parents=True, exist_ok=True)
    source = ROOT / "inputs" / "ligands" / f"{item['ligand_id']}_PubChem3D.sdf"
    protonated = outdir / f"{item['ligand_id']}_pH7.4.sdf"
    ligand_pdbqt = outdir / f"{item['ligand_id']}.pdbqt"
    docked_pdbqt = outdir / f"{item['ligand_id']}_docked.pdbqt"
    docked_sdf = outdir / f"{item['ligand_id']}_docked.sdf"

    run(["obabel", "-isdf", str(source), "-osdf", "-O", str(protonated), "-p", "7.4"],
        LOGS / f"03_{target}_{item['ligand_id']}_protonation.log")
    run(["mk_prepare_ligand.py", "-i", str(protonated), "-o", str(ligand_pdbqt),
         "--charge_model", "gasteiger", "--add_index_map"],
        LOGS / f"03_{target}_{item['ligand_id']}_ligand_prep.log")
    run([
        "vina", "--receptor", str(PREP / target / "receptor.pdbqt"),
        "--ligand", str(ligand_pdbqt), "--out", str(docked_pdbqt),
        "--center_x", box["center_x"], "--center_y", box["center_y"], "--center_z", box["center_z"],
        "--size_x", box["size_x"], "--size_y", box["size_y"], "--size_z", box["size_z"],
        "--exhaustiveness", str(EXHAUSTIVENESS), "--seed", str(SEED),
        "--num_modes", str(NUM_MODES), "--energy_range", "8", "--cpu", "8",
    ], LOGS / f"03_{target}_{item['ligand_id']}_vina.log")
    run(["mk_export.py", str(docked_pdbqt), "-s", str(docked_sdf)],
        LOGS / f"03_{target}_{item['ligand_id']}_export.log")

    scores = affinities(docked_pdbqt)
    if not scores:
        raise RuntimeError(f"No Vina score parsed for {target}/{item['ligand_id']}")
    rows.append({
        **item,
        "ligand_source": f"PubChem CID {item['pubchem_cid']} 3D conformer; identity cross-checked to TCMSP {item['ligand_id']}",
        "target_ligand_mapping_source": "project archive 02.tcmsp_target_dangshen.xlsx, Sheet2",
        "pdb_id": box["pdb_id"],
        "chain": box["chain"],
        "binding_site_definition": f"co-crystallized {box['native_ligand']} pocket",
        "grid_center": f"{box['center_x']};{box['center_y']};{box['center_z']}",
        "grid_size_A": f"{box['size_x']};{box['size_y']};{box['size_z']}",
        "exhaustiveness": EXHAUSTIVENESS,
        "random_seed": SEED,
        "num_modes_generated": len(scores),
        "best_affinity_kcal_mol": scores[0],
        "interpretation": "in-silico pose compatibility only; not evidence of physical binding or efficacy",
        "best_pose_file": str(docked_sdf.relative_to(ROOT)),
    })

with (ROOT / "new_docking_results.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
