#!/usr/bin/env python3
"""Report geometric contacts for the final revised docking results."""

from pathlib import Path
import csv
import math

from rdkit import Chem


ROOT = Path(__file__).resolve().parents[1]
POLAR = {"N", "O", "S"}


def distance(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def receptor_atoms(path):
    atoms = []
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        element = line[76:78].strip().upper()
        if not element or element == "H":
            continue
        atoms.append({
            "atom": line[12:16].strip(), "resname": line[17:20].strip(),
            "chain": line[21].strip(), "resnum": int(line[22:26]), "element": element,
            "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
        })
    return atoms


with (ROOT / "new_docking_results.csv").open() as handle:
    rows = list(csv.DictReader(handle))

interaction_rows = []
for row in rows:
    receptor = receptor_atoms(ROOT / "prepared" / row["target"] / "receptor_pH7.4.pdb")
    pose_path = Path(row["best_pose_file"])
    if not pose_path.is_absolute():
        pose_path = ROOT / pose_path
    pose = next(mol for mol in Chem.SDMolSupplier(str(pose_path), removeHs=False)
                if mol is not None)
    conf = pose.GetConformer()
    ligand_atoms = []
    for atom in pose.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        point = conf.GetAtomPosition(atom.GetIdx())
        ligand_atoms.append({"idx": atom.GetIdx() + 1, "element": atom.GetSymbol().upper(),
                             "xyz": (point.x, point.y, point.z)})

    residue_min, polar_pairs = {}, []
    for ligand_atom in ligand_atoms:
        for receptor_atom in receptor:
            d = distance(ligand_atom["xyz"], receptor_atom["xyz"])
            if receptor_atom["resname"] in {"ZN", "CA"}:
                continue
            key = f"{receptor_atom['resname']}{receptor_atom['resnum']}({receptor_atom['chain']})"
            residue_min[key] = min(residue_min.get(key, 999.0), d)
            if (ligand_atom["element"] in POLAR and receptor_atom["element"] in POLAR and d <= 3.5):
                polar_pairs.append(
                    f"L{ligand_atom['element']}{ligand_atom['idx']}-{key}:{receptor_atom['atom']} {d:.2f} A"
                )
    contacts = sorted(((key, value) for key, value in residue_min.items() if value <= 4.0),
                      key=lambda item: item[1])
    zinc = [atom for atom in receptor if atom["resname"] == "ZN"]
    hetero = [atom for atom in ligand_atoms if atom["element"] in POLAR]
    zinc_distance = min((distance(a["xyz"], z["xyz"]) for a in hetero for z in zinc), default=None)
    interaction_rows.append({
        "target": row["target"], "ligand_id": row["ligand_id"], "ligand": row["ligand"],
        "contact_residues_within_4A": "; ".join(f"{key}:{value:.2f} A" for key, value in contacts),
        "polar_heavy_atom_contacts_within_3.5A": "; ".join(sorted(set(polar_pairs))),
        "minimum_ligand_heteroatom_to_Zn_distance_A": "" if zinc_distance is None else f"{zinc_distance:.3f}",
        "interpretation": "distance-based geometric contacts; polar proximity is not asserted as a hydrogen bond",
    })

with (ROOT / "new_docking_interactions.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=interaction_rows[0].keys())
    writer.writeheader()
    writer.writerows(interaction_rows)
