#!/usr/bin/env python3
"""Summarize reproducible geometric contacts for each top-ranked pose."""

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
    results = list(csv.DictReader(handle))

interaction_rows = []
for row in results:
    target = row["target"]
    rec_atoms = receptor_atoms(ROOT / "prepared" / target / "receptor_pH7.4.pdb")
    mol = next(m for m in Chem.SDMolSupplier(str(ROOT / row["best_pose_file"]), removeHs=False) if m is not None)
    conf = mol.GetConformer()
    lig_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        p = conf.GetAtomPosition(atom.GetIdx())
        lig_atoms.append({"idx": atom.GetIdx(), "element": atom.GetSymbol().upper(), "xyz": (p.x, p.y, p.z)})

    residue_min = {}
    polar_pairs = []
    for la in lig_atoms:
        for ra in rec_atoms:
            d = distance(la["xyz"], ra["xyz"])
            if ra["resname"] != "ZN":
                key = f"{ra['resname']}{ra['resnum']}({ra['chain']})"
                residue_min[key] = min(residue_min.get(key, 999.0), d)
                if la["element"] in POLAR and ra["element"] in POLAR and d <= 3.5:
                    polar_pairs.append(f"L{la['element']}{la['idx'] + 1}-{key}:{ra['atom']} {d:.2f}A")

    contacts = sorted(((k, v) for k, v in residue_min.items() if v <= 4.0), key=lambda x: x[1])
    zn_atoms = [a for a in rec_atoms if a["resname"] == "ZN"]
    zn_dist = ""
    if zn_atoms:
        donors = [a for a in lig_atoms if a["element"] in POLAR]
        if donors:
            zn_dist = min(distance(a["xyz"], z["xyz"]) for a in donors for z in zn_atoms)

    contact_text = "; ".join(f"{k}:{v:.2f}A" for k, v in contacts)
    polar_text = "; ".join(sorted(set(polar_pairs)))
    row["contact_residues_within_4A"] = contact_text
    row["polar_heavy_atom_contacts_within_3.5A"] = polar_text
    row["minimum_ligand_heteroatom_to_Zn_distance_A"] = "" if zn_dist == "" else f"{zn_dist:.3f}"
    row["interaction_definition"] = "distance-based geometric contacts; polar contacts are not asserted as hydrogen bonds"
    interaction_rows.append({
        "target": target, "ligand_id": row["ligand_id"], "ligand": row["ligand"],
        "contact_residues_within_4A": contact_text,
        "polar_heavy_atom_contacts_within_3.5A": polar_text,
        "minimum_ligand_heteroatom_to_Zn_distance_A": row["minimum_ligand_heteroatom_to_Zn_distance_A"],
    })

with (ROOT / "new_docking_results.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=results[0].keys())
    writer.writeheader()
    writer.writerows(results)
with (ROOT / "new_docking_interactions.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=interaction_rows[0].keys())
    writer.writeheader()
    writer.writerows(interaction_rows)
