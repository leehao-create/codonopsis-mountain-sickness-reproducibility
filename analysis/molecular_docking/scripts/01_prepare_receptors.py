#!/usr/bin/env python3
"""Extract receptors/native ligands and prepare reproducible Vina inputs."""

from pathlib import Path
import csv
import math
import subprocess
import shutil

from rdkit import Chem
from rdkit.Chem import AllChem
from pdbfixer import PDBFixer
from openmm.app import PDBFile


ROOT = Path(__file__).resolve().parents[1]
PREP = ROOT / "prepared"
LOGS = ROOT / "logs"
PREP.mkdir(parents=True, exist_ok=True)
LOGS.mkdir(exist_ok=True)

TARGETS = {
    "NR3C2": {
        "pdb_id": "2AA2", "chain": "A", "resolution_A": 1.95,
        "native_resname": "AS4", "native_resnum": 201,
        "native_name": "aldosterone", "component": "AS4_ideal.sdf",
        "retain_hetero": set(),
        "construct": "human NR3C2 ligand-binding domain; engineered C808S mutation; no S810L",
        "wild_type": "no (C808S crystallization mutation)",
    },
    "MMP9": {
        "pdb_id": "1GKC", "chain": "A", "resolution_A": 2.30,
        "native_resname": "NFH", "native_resnum": 1448,
        "native_name": "reverse-hydroxamate inhibitor NFH", "component": "NFH_ideal.sdf",
        "retain_hetero": {("CA", 1444), ("CA", 1445), ("CA", 1446), ("CA", 1447),
                            ("CA", 1452), ("ZN", 1450), ("ZN", 1451)},
        "construct": "wild-type human MMP9 catalytic domain, residues 107-215 and 391-443",
        "wild_type": "yes (no engineered point mutation reported)",
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
    return [value[1] for value in chosen.values()]


def run(cmd, log_path):
    result = subprocess.run(cmd, text=True, capture_output=True)
    log_path.write_text("COMMAND: " + " ".join(map(str, cmd)) + "\n\n" + result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f"Command failed: {' '.join(map(str, cmd))}\n{result.stderr}")


rows = []
for target, cfg in TARGETS.items():
    target_dir = PREP / target
    target_dir.mkdir(exist_ok=True)
    pdb_file = ROOT / "inputs" / "structures" / f"{cfg['pdb_id']}.pdb"
    all_lines = pdb_file.read_text().splitlines()

    receptor_lines, ligand_lines, ligand_serials = [], [], set()
    for line in all_lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        _, resname, chain, resnum = atom_identity(line)
        if chain != cfg["chain"]:
            continue
        if line.startswith("ATOM  "):
            receptor_lines.append(line)
        elif resname == cfg["native_resname"] and resnum == cfg["native_resnum"]:
            ligand_lines.append(line)
            ligand_serials.add(int(line[6:11]))
        elif (resname, resnum) in cfg["retain_hetero"]:
            receptor_lines.append(line)

    receptor_lines = choose_altloc(receptor_lines)
    ligand_lines = choose_altloc(ligand_lines)
    if not ligand_lines:
        raise RuntimeError(f"Native ligand not found for {target}")

    receptor_raw = target_dir / "receptor_raw_chainA_metals_retained.pdb"
    receptor_raw.write_text("\n".join(receptor_lines + ["TER", "END"]) + "\n")

    fixer = PDBFixer(filename=str(receptor_raw))
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    receptor_h = target_dir / "receptor_pH7.4.pdb"
    with receptor_h.open("w") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)

    # Open Babel is used for the rigid protein PDBQT because it retains the
    # complete PDBFixer-prepared C terminus. MMP9 metals are appended below
    # with their original coordinates, +2 charges and Vina Ca/Zn atom types.
    receptor_for_obabel = target_dir / "receptor_for_openbabel.pdb"
    receptor_h_lines = receptor_h.read_text().splitlines()
    metal_lines = [line for line in receptor_h_lines
                   if line.startswith("HETATM") and line[17:20].strip() in {"ZN", "CA"}]
    receptor_for_obabel.write_text("\n".join(
        line for line in receptor_h_lines
        if not (line.startswith("HETATM") and line[17:20].strip() in {"ZN", "CA"})
    ) + "\n")
    receptor_pdbqt = target_dir / "receptor.pdbqt"
    run([
        shutil.which("obabel") or "obabel", "-ipdb", str(receptor_for_obabel), "-opdbqt",
        "-O", str(receptor_pdbqt), "-xr",
        "--partialcharge", "gasteiger",
    ], LOGS / f"01_{target}_receptor_openbabel.log")
    if receptor_pdbqt.stat().st_size == 0:
        raise RuntimeError(f"Open Babel generated an empty receptor PDBQT for {target}")
    if metal_lines:
        pdbqt_lines = [line for line in receptor_pdbqt.read_text().splitlines()
                       if not line.startswith(("TER", "END"))]
        serial = sum(line.startswith(("ATOM", "HETATM")) for line in pdbqt_lines)
        for line in metal_lines:
            serial += 1
            atom = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21].strip()
            resnum = int(line[22:26])
            x, y, z = xyz(line)
            atom_type = "Zn" if resname == "ZN" else "Ca"
            pdbqt_lines.append(
                f"HETATM{serial:5d} {atom:<4s} {resname:>3s} {chain}{resnum:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00    +2.000 {atom_type:>2s}"
            )
        receptor_pdbqt.write_text("\n".join(pdbqt_lines + ["TER"]) + "\n")

    conect = []
    for line in all_lines:
        if not line.startswith("CONECT"):
            continue
        serials = [int(line[i:i + 5]) for i in range(6, len(line), 5) if line[i:i + 5].strip()]
        if serials and serials[0] in ligand_serials:
            filtered = [serial for serial in serials if serial in ligand_serials]
            if len(filtered) > 1:
                conect.append("CONECT" + "".join(f"{serial:5d}" for serial in filtered))
    ligand_pdb = target_dir / "native_ligand_crystal.pdb"
    ligand_pdb.write_text("\n".join(ligand_lines + conect + ["END"]) + "\n")

    template = Chem.SDMolSupplier(
        str(ROOT / "inputs" / "native_ligands" / cfg["component"]), removeHs=False
    )[0]
    pdb_mol = Chem.MolFromPDBFile(str(ligand_pdb), removeHs=False, sanitize=False, proximityBonding=True)
    if template is None or pdb_mol is None:
        raise RuntimeError(f"Cannot parse native ligand for {target}")
    crystal = AllChem.AssignBondOrdersFromTemplate(Chem.RemoveHs(template), Chem.RemoveHs(pdb_mol))
    Chem.SanitizeMol(crystal)
    crystal.SetProp("_Name", f"{cfg['pdb_id']}_{cfg['native_resname']}_crystal")
    writer = Chem.SDWriter(str(target_dir / "native_ligand_crystal.sdf"))
    writer.write(crystal)
    writer.close()

    coords = [xyz(line) for line in ligand_lines]
    center = [sum(point[i] for point in coords) / len(coords) for i in range(3)]
    extent = [max(point[i] for point in coords) - min(point[i] for point in coords) for i in range(3)]
    box = [max(22.0, float(math.ceil(value + 12.0))) for value in extent]
    rows.append({
        "target": target, "pdb_id": cfg["pdb_id"], "chain": cfg["chain"],
        "resolution_A": cfg["resolution_A"], "construct": cfg["construct"],
        "wild_type": cfg["wild_type"], "co_crystallized_ligand": cfg["native_name"],
        "co_crystallized_ligand_id": cfg["native_resname"],
        "docking_site_definition": "co-crystallized ligand heavy-atom centroid",
        "center_x": round(center[0], 3), "center_y": round(center[1], 3),
        "center_z": round(center[2], 3), "size_x": box[0], "size_y": box[1],
        "size_z": box[2],
        "box_rule": "native-ligand heavy-atom extent plus 12 A total padding; minimum 22 A per axis",
        "retained_nonprotein_components": "none" if not cfg["retain_hetero"] else "all chain-A Zn2+ and Ca2+ ions",
        "receptor_pdbqt_method": (
            "Open Babel 3.2.1 rigid protein with Gasteiger partial charges; "
            "original-coordinate metals appended as +2 Ca/Zn Vina types" if metal_lines else
            "Open Babel 3.2.1 rigid receptor with Gasteiger partial charges"
        ),
    })

with (PREP / "receptor_and_grid_metadata.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
