#!/usr/bin/env python3
"""Extract one receptor chain and its crystallographic ligand reproducibly."""

from pathlib import Path
import csv
import math

from rdkit import Chem
from rdkit.Chem import AllChem
from pdbfixer import PDBFixer
from openmm.app import PDBFile


ROOT = Path(__file__).resolve().parents[1]
STRUCTURES = ROOT / "inputs" / "structures"
LIGANDS = ROOT / "inputs" / "ligands"
PREP = ROOT / "prepared"
PREP.mkdir(parents=True, exist_ok=True)

TARGETS = {
    "CA2": {
        "pdb_id": "5E28",
        "chain": "A",
        "native_resname": "BC5",
        "native_resnum": 302,
        "component": "BC5_ideal.sdf",
        "retain_hetero": {("ZN", 301)},
        "resolution": 1.30,
        "construct": "human CA2 residues 4-260; no engineered mutation",
    },
    "BCL2": {
        "pdb_id": "8HTS",
        "chain": "A",
        "native_resname": "N2L",
        "native_resnum": 301,
        "component": "N2L_ideal.sdf",
        "retain_hetero": set(),
        "resolution": 1.25,
        "construct": "engineered human BCL2 construct; RCSB reports no engineered point mutation",
    },
}


def atom_identity(line):
    return line[12:16].strip(), line[17:20].strip(), line[21].strip(), int(line[22:26])


def xyz(line):
    return tuple(float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54)))


def choose_altloc(lines):
    chosen = {}
    for line in lines:
        key = atom_identity(line)
        alt = line[16].strip()
        occ = float(line[54:60] or 0)
        rank = (alt not in ("", "A"), -occ, alt)
        if key not in chosen or rank < chosen[key][0]:
            chosen[key] = (rank, line[:16] + " " + line[17:])
    return [v[1] for v in chosen.values()]


rows = []
for target, cfg in TARGETS.items():
    target_dir = PREP / target
    target_dir.mkdir(exist_ok=True)
    pdb_file = STRUCTURES / f"{cfg['pdb_id']}.pdb"
    all_lines = pdb_file.read_text().splitlines()

    atom_lines = []
    ligand_lines = []
    ligand_serials = set()
    for line in all_lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        _, resname, chain, resnum = atom_identity(line)
        if chain != cfg["chain"]:
            continue
        if line.startswith("ATOM  "):
            atom_lines.append(line)
        elif resname == cfg["native_resname"] and resnum == cfg["native_resnum"]:
            ligand_lines.append(line)
            ligand_serials.add(int(line[6:11]))
        elif (resname, resnum) in cfg["retain_hetero"]:
            atom_lines.append(line)

    atom_lines = choose_altloc(atom_lines)
    ligand_lines = choose_altloc(ligand_lines)
    if not ligand_lines:
        raise RuntimeError(f"No native ligand found for {target}")

    receptor_file = target_dir / "receptor_raw.pdb"
    receptor_file.write_text("\n".join(atom_lines + ["TER", "END"]) + "\n")

    # Repair only missing heavy atoms and add hydrogens at pH 7.4. Missing
    # residues are intentionally not modeled into an experimental construct.
    fixer = PDBFixer(filename=str(receptor_file))
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    with (target_dir / "receptor_pH7.4.pdb").open("w") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)

    conect = []
    for line in all_lines:
        if line.startswith("CONECT"):
            serials = [int(line[i:i + 5]) for i in range(6, len(line), 5) if line[i:i + 5].strip()]
            if serials and serials[0] in ligand_serials:
                filtered = [s for s in serials if s in ligand_serials]
                if len(filtered) > 1:
                    conect.append("CONECT" + "".join(f"{s:5d}" for s in filtered))
    ligand_pdb = target_dir / "native_ligand_crystal.pdb"
    ligand_pdb.write_text("\n".join(ligand_lines + conect + ["END"]) + "\n")

    template = Chem.SDMolSupplier(str(LIGANDS / cfg["component"]), removeHs=False)[0]
    pdb_mol = Chem.MolFromPDBFile(str(ligand_pdb), removeHs=False, sanitize=False, proximityBonding=True)
    if template is None or pdb_mol is None:
        raise RuntimeError(f"Failed reading native ligand for {target}")
    try:
        crystal = AllChem.AssignBondOrdersFromTemplate(Chem.RemoveHs(template), Chem.RemoveHs(pdb_mol))
    except Exception as exc:
        raise RuntimeError(f"Bond-order assignment failed for {target}: {exc}") from exc
    Chem.SanitizeMol(crystal)
    crystal.SetProp("_Name", f"{cfg['pdb_id']}_{cfg['native_resname']}_crystal")
    writer = Chem.SDWriter(str(target_dir / "native_ligand_crystal.sdf"))
    writer.write(crystal)
    writer.close()

    coords = [xyz(line) for line in ligand_lines]
    center = [sum(v[i] for v in coords) / len(coords) for i in range(3)]
    extent = [max(v[i] for v in coords) - min(v[i] for v in coords) for i in range(3)]
    # Predefined rule: ligand extent plus 12 A total padding, minimum 22 A.
    box = [max(22.0, math.ceil(e + 12.0)) for e in extent]

    rows.append({
        "target": target,
        "pdb_id": cfg["pdb_id"],
        "chain": cfg["chain"],
        "resolution_A": cfg["resolution"],
        "construct": cfg["construct"],
        "native_ligand": cfg["native_resname"],
        "center_x": round(center[0], 3),
        "center_y": round(center[1], 3),
        "center_z": round(center[2], 3),
        "size_x": box[0],
        "size_y": box[1],
        "size_z": box[2],
        "box_rule": "native-ligand heavy-atom extent + 12 A total padding; minimum 22 A",
    })

with (PREP / "docking_boxes.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
