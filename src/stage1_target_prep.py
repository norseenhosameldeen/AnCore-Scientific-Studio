#!/usr/bin/env python3
"""
AnCore-Scientific-Studio v1.2 — Stage 1: Target Receptor Curation,
PDBFixer + OpenMM AMBER14 Sidechain Relaxation, and 5.0 Å Pocket Box Mapping.

Save location in your local project:
    AnCore_v1.2/src/stage1_target_prep.py

CLI usage:
  1) Run all 5 benchmark targets (1IEP, 3ERT, 1UYG, 1M17, 1STP) in one command (DEFAULT):
     python3 src/stage1_target_prep.py --run-all

  2) Include the de novo custom target 4WKQ as well (all 6 paper targets):
     python3 src/stage1_target_prep.py --run-all --include-4wkq

  3) Run a single custom or benchmark target:
     python3 src/stage1_target_prep.py --pdb-id 4WKQ --out-dir workspaces/4WKQ
"""

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# ── Built-In Multi-Target Benchmark Suite (Diverse Conformational & Family States) ──
BENCHMARK_SUITE = [
    {
        "pdb_id": "1IEP",
        "target_name": "ABL1 tyrosine kinase (DFG-out inactive tunnel)",
        "preferred_chain": "A",
        "preferred_ligand": "STI",
    },
    {
        "pdb_id": "3ERT",
        "target_name": "Estrogen receptor alpha (antagonist nuclear receptor LBD)",
        "preferred_chain": "A",
        "preferred_ligand": "OHT",
    },
    {
        "pdb_id": "1UYG",
        "target_name": "HSP90-alpha ATPase domain (polar Bergerat fold)",
        "preferred_chain": "A",
        "preferred_ligand": "PU2",
    },
    {
        "pdb_id": "1M17",
        "target_name": "EGFR tyrosine kinase (DFG-in active hinge cleft)",
        "preferred_chain": "A",
        "preferred_ligand": "AQ4",
    },
    {
        "pdb_id": "1STP",
        "target_name": "Streptavidin beta-barrel (femtomolar polar H-bond cage)",
        "preferred_chain": "A",
        "preferred_ligand": "BTN",
    },
]

CUSTOM_4WKQ_TARGET = {
    "pdb_id": "4WKQ",
    "target_name": "EGFR kinase domain bound to gefitinib (de novo custom case)",
    "preferred_chain": "A",
    "preferred_ligand": "IRE",
}

# ── Standard Amino Acid & Exclusion Tables ────────────────────────────────────
STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP", "CYX", "ASH", "GLH", "LYN",
}

EXCLUDED_HETATM = {
    "HOH", "WAT", "TIP", "SOL", "DOD", "SO4", "PO4", "GOL", "EDO", "DMS",
    "ACT", "ACE", "NME", "NAG", "NDG", "BMA", "MAN", "GAL", "FUC", "PEG",
    "PGE", "PG4", "MPD", "TRS", "MES", "EPE", "HEPES", "CIT", "TAR", "FMT",
    "NA", "CL", "K", "MG", "CA", "ZN", "FE", "MN", "CU", "CO", "NI", "CD",
    "IOD", "BR", "HG", "URE", "SUL", "NO3", "CO3", "BCT", "UNX", "UNL",
}

RESIDUE_CLASS_MAP = {
    "PHE": "Aromatic", "TYR": "Aromatic", "TRP": "Aromatic", "HIS": "Aromatic",
    "HID": "Aromatic", "HIE": "Aromatic", "HIP": "Aromatic",
    "ASP": "Charged/Salt-Bridge", "GLU": "Charged/Salt-Bridge",
    "ARG": "Charged/Salt-Bridge", "LYS": "Charged/Salt-Bridge",
    "SER": "Polar H-Bonding", "THR": "Polar H-Bonding",
    "ASN": "Polar H-Bonding", "GLN": "Polar H-Bonding", "CYS": "Polar H-Bonding",
    "ALA": "Hydrophobic", "VAL": "Hydrophobic", "LEU": "Hydrophobic",
    "ILE": "Hydrophobic", "MET": "Hydrophobic", "PRO": "Hydrophobic",
    "GLY": "Hydrophobic",
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def format_strict_pdb_hetatm(
    atom_serial: int,
    atom_name: str,
    x: float,
    y: float,
    z: float,
    element: str,
    res_name: str = "LIG",
    chain_id: str = "L",
    res_seq: int = 1,
) -> str:
    """Writes a strictly compliant 80-column PDB HETATM record (Fixes Trap 4)."""
    el = element.strip().upper()[:2] or "C"
    name_clean = atom_name.strip()[:4]
    name_field = f" {name_clean:<3s}" if len(el) == 1 and len(name_clean) < 4 else f"{name_clean:<4s}"
    return (
        f"HETATM{atom_serial:5d} {name_field:4s} {res_name:>3s} {chain_id:1s}{res_seq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}  \n"
    )


def fetch_or_locate_raw_pdb(
    pdb_id: str,
    local_pdb: Optional[Path],
    out_dir: Path,
    search_dirs: Optional[List[Path]] = None,
) -> Path:
    """Locates a local .pdb file or downloads it live from RCSB PDB."""
    out_dir.mkdir(parents=True, exist_ok=True)
    target_raw = out_dir / f"{pdb_id.upper()}_raw.pdb"

    if local_pdb is not None and local_pdb.is_file():
        if local_pdb.resolve() != target_raw.resolve():
            target_raw.write_bytes(local_pdb.read_bytes())
        return target_raw

    candidates = [target_raw]
    if search_dirs:
        for sdir in search_dirs:
            candidates.extend([
                sdir / f"{pdb_id.upper()}.pdb",
                sdir / f"{pdb_id.lower()}.pdb",
                sdir / pdb_id.upper() / f"{pdb_id.upper()}_raw.pdb",
                sdir / pdb_id.upper() / "receptor_raw.pdb",
                sdir / pdb_id.upper() / f"{pdb_id.upper()}.pdb",
            ])
    for cand in candidates:
        if cand.is_file() and cand.stat().st_size > 100:
            if cand.resolve() != target_raw.resolve():
                target_raw.write_bytes(cand.read_bytes())
            return target_raw

    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    print(f"  [RCSB Fetch] Downloading {pdb_id.upper()} from {url} ...")
    req = urllib.request.Request(url, headers={"User-Agent": "AnCore-Scientific-Studio/1.2"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        data = resp.read()
    target_raw.write_bytes(data)
    return target_raw


def isolate_primary_chain_and_ligand(
    raw_pdb: Path,
    preferred_chain: Optional[str] = None,
    preferred_ligand: Optional[str] = None,
) -> Tuple[str, str, List[str], List[dict], List[dict]]:
    """Isolates the primary receptor chain and primary co-crystallized ligand (heavy atoms)."""
    lines = raw_pdb.read_text(encoding="utf-8", errors="ignore").splitlines()

    chain_atoms: Dict[str, List[dict]] = {}
    chain_lines: Dict[str, List[str]] = {}
    het_groups: Dict[Tuple[str, str, int], List[dict]] = {}

    for line in lines:
        rec = line[:6].strip()
        if rec not in ("ATOM", "HETATM") or len(line) < 54:
            continue
        alt_loc = line[16:17].strip()
        if alt_loc not in ("", "A", "1"):
            continue

        atom_name = line[12:16].strip()
        res_name = line[17:20].strip().upper()
        chain_id = line[21:22].strip() or "A"
        try:
            res_seq = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
        except ValueError:
            continue

        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if not element:
            stripped = "".join(c for c in atom_name if c.isalpha())
            element = stripped[:1].upper() if stripped else "C"

        if element in ("H", "D") or atom_name.upper().startswith(("H", "1H", "2H", "3H")):
            continue

        atom_info = {
            "line": line if line.endswith("\n") else line + "\n",
            "atom_name": atom_name,
            "res_name": res_name,
            "chain_id": chain_id,
            "res_seq": res_seq,
            "x": x,
            "y": y,
            "z": z,
            "element": element,
        }

        if rec == "ATOM" and res_name in STANDARD_AA:
            chain_atoms.setdefault(chain_id, []).append(atom_info)
            chain_lines.setdefault(chain_id, []).append(atom_info["line"])
        elif rec == "HETATM" and res_name not in EXCLUDED_HETATM and res_name not in STANDARD_AA:
            key = (res_name, chain_id, res_seq)
            het_groups.setdefault(key, []).append(atom_info)

    if not chain_atoms:
        raise RuntimeError(f"No standard protein ATOM records found in {raw_pdb}")

    if preferred_chain and preferred_chain in chain_atoms:
        sel_chain = preferred_chain
    else:
        sel_chain = max(chain_atoms.keys(), key=lambda c: len(chain_atoms[c]))

    rec_atoms = chain_atoms[sel_chain]
    rec_coords = np.array([[a["x"], a["y"], a["z"]] for a in rec_atoms])
    rec_centroid = rec_coords.mean(axis=0)

    best_lig_key = None
    best_score = -1e9
    for key, l_atoms in het_groups.items():
        res_name, l_chain, _ = key
        n_heavy = len(l_atoms)
        if n_heavy < 5:
            continue
        if preferred_ligand and preferred_ligand.upper() != "AUTO":
            if res_name != preferred_ligand.upper():
                continue
        l_coords = np.array([[a["x"], a["y"], a["z"]] for a in l_atoms])
        min_d = float(np.linalg.norm(rec_coords[:, None, :] - l_coords[None, :, :], axis=-1).min())
        if min_d > 6.5:
            continue
        chain_bonus = 50.0 if l_chain == sel_chain else 0.0
        score = n_heavy + chain_bonus - min_d
        if score > best_score:
            best_score = score
            best_lig_key = key

    if best_lig_key is not None:
        lig_res_name = best_lig_key[0]
        lig_atoms = het_groups[best_lig_key]
    else:
        lig_res_name = "LIG"
        lig_atoms = [{
            "atom_name": "C1",
            "res_name": "LIG",
            "chain_id": sel_chain,
            "res_seq": 1,
            "x": float(rec_centroid[0]),
            "y": float(rec_centroid[1]),
            "z": float(rec_centroid[2]),
            "element": "C",
        }]

    return sel_chain, lig_res_name, chain_lines[sel_chain], rec_atoms, lig_atoms


def repair_and_minimize_receptor(
    raw_chain_pdb: Path,
    minimized_pdb: Path,
    ph: float = 7.4,
    max_iters: int = 200,
) -> dict:
    """Runs PDBFixer + 200-step OpenMM AMBER14 sidechain minimization with backbone restraints."""
    info = {
        "pdbfixer_applied": False,
        "openmm_minimized": False,
        "forcefield": "amber14-all.xml + implicit/gbn2.xml",
        "ph": ph,
        "max_iterations": max_iters,
    }
    try:
        from pdbfixer import PDBFixer
        import openmm as mm
        import openmm.app as app
        import openmm.unit as unit

        fixer = PDBFixer(filename=str(raw_chain_pdb))
        fixer.findMissingResidues()
        chains = list(fixer.topology.chains())
        keys_to_remove = []
        for key, res_list in fixer.missingResidues.items():
            chain_idx, res_idx = key
            if res_idx == 0 or (chain_idx < len(chains) and res_idx >= len(list(chains[chain_idx].residues()))) or len(res_list) > 4:
                keys_to_remove.append(key)
        for k in keys_to_remove:
            del fixer.missingResidues[k]

        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(keepWater=False)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(ph)
        info["pdbfixer_applied"] = True

        forcefield = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
        system = forcefield.createSystem(
            fixer.topology,
            nonbondedMethod=app.CutoffNonPeriodic,
            nonbondedCutoff=1.2 * unit.nanometer,
            constraints=app.HBonds,
        )

        restraint = mm.CustomExternalForce("0.5 * k_ca * ((x - x0)^2 + (y - y0)^2 + (z - z0)^2)")
        restraint.addGlobalParameter("k_ca", 10.0 * 418.4)
        restraint.addPerParticleParameter("x0")
        restraint.addPerParticleParameter("y0")
        restraint.addPerParticleParameter("z0")

        pos_nm = fixer.positions.value_in_unit(unit.nanometer)
        for atom in fixer.topology.atoms():
            if atom.name in ("CA", "N", "C"):
                p = pos_nm[atom.index]
                restraint.addParticle(atom.index, [p.x, p.y, p.z])
        system.addForce(restraint)

        integrator = mm.LangevinMiddleIntegrator(300.0 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds)
        simulation = app.Simulation(fixer.topology, system, integrator)
        simulation.context.setPositions(fixer.positions)

        state_pre = simulation.context.getState(getEnergy=True)
        e_pre = state_pre.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)

        simulation.minimizeEnergy(maxIterations=max_iters)

        state_post = simulation.context.getState(getPositions=True, getEnergy=True)
        e_post = state_post.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)

        with open(minimized_pdb, "w", encoding="utf-8") as f:
            app.PDBFile.writeFile(fixer.topology, state_post.getPositions(), f, keepIds=True)

        info["openmm_minimized"] = True
        info["energy_pre_kcal_mol"] = round(float(e_pre), 2)
        info["energy_post_kcal_mol"] = round(float(e_post), 2)
        return info

    except Exception as exc:
        raise RuntimeError(
            f"[Stage 1 Physics Error] PDBFixer / OpenMM sidechain minimization failed: {exc}."
        ) from exc


def prepare_receptor_pdbqt(minimized_pdb: Path, out_pdbqt: Path) -> str:
    """Converts the pH 7.4 minimized receptor PDB into a rigid AutoDock Vina .pdbqt file."""
    obabel_bin = shutil.which("obabel")
    if obabel_bin:
        cmd = [obabel_bin, "-ipdb", str(minimized_pdb), "-opdbqt", "-O", str(out_pdbqt), "-xr"]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode == 0 and out_pdbqt.is_file() and out_pdbqt.stat().st_size > 500:
            return "OpenBabel (-xr rigid receptor)"

    AROMATIC_ATOMS = {
        ("PHE", "CG"), ("PHE", "CD1"), ("PHE", "CD2"), ("PHE", "CE1"), ("PHE", "CE2"), ("PHE", "CZ"),
        ("TYR", "CG"), ("TYR", "CD1"), ("TYR", "CD2"), ("TYR", "CE1"), ("TYR", "CE2"), ("TYR", "CZ"),
        ("TRP", "CG"), ("TRP", "CD1"), ("TRP", "CD2"), ("TRP", "CE2"), ("TRP", "CE3"), ("TRP", "CZ2"), ("TRP", "CZ3"), ("TRP", "CH2"),
        ("HIS", "CG"), ("HIS", "CD2"), ("HIS", "CE1"),
        ("HID", "CG"), ("HID", "CD2"), ("HID", "CE1"),
        ("HIE", "CG"), ("HIE", "CD2"), ("HIE", "CE1"),
        ("HIP", "CG"), ("HIP", "CD2"), ("HIP", "CE1"),
    }

    lines_out = []
    for line in minimized_pdb.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith("ATOM  "):
            continue
        atom_name = line[12:16].strip()
        res_name = line[17:20].strip().upper()
        el = line[76:78].strip().upper() if len(line) >= 78 else ""
        if not el:
            el = "".join(c for c in atom_name if c.isalpha())[:1].upper() or "C"

        if el == "H":
            if not any(atom_name.startswith(p) for p in ("H", "HN", "HG", "HH", "HD", "HE", "HZ")):
                continue
            if atom_name.startswith(("HA", "HB", "HG1", "HG2", "HD1", "HD2", "HE1", "HE2", "HE3")):
                if res_name not in ("HIS", "HID", "HIE", "HIP", "TRP", "TYR", "SER", "THR", "CYS", "ASN", "GLN", "LYS", "ARG"):
                    continue
            ad_type = "HD"
            charge = 0.160
        elif el == "C":
            ad_type = "A" if (res_name, atom_name) in AROMATIC_ATOMS else "C"
            charge = 0.050
        elif el == "O":
            ad_type = "OA"
            charge = -0.550 if atom_name in ("O", "OD1", "OD2", "OE1", "OE2") else -0.390
        elif el == "N":
            ad_type = "NA" if (res_name in ("HIS", "HID", "HIE") and atom_name in ("ND1", "NE2")) else "N"
            charge = -0.350
        elif el == "S":
            ad_type = "SA" if res_name == "CYS" else "S"
            charge = -0.180
        else:
            ad_type = el[:2].capitalize()
            charge = 0.000

        base_cols = line[:54]
        lines_out.append(f"{base_cols}  1.00  0.00    {charge:+6.3f} {ad_type:<2s}\n")

    out_pdbqt.write_text("".join(lines_out), encoding="utf-8")
    return "AnCore AD4 United-Atom Receptor Compiler"


def compute_pocket_and_box(
    raw_rec_atoms: List[dict],
    min_rec_pdb: Path,
    lig_atoms: List[dict],
    pocket_cutoff_A: float = 5.0,
    box_padding_A: float = 4.0,
    s_min_A: float = 16.0,
    s_max_A: float = 24.0,
) -> dict:
    """Computes the 5.0 Å pocket (Eq. 1) and 3D search box (Eqs. 2–3)."""
    lig_coords = np.array([[a["x"], a["y"], a["z"]] for a in lig_atoms], dtype=float)

    res_groups: Dict[Tuple[str, int, str], List[np.ndarray]] = {}
    for a in raw_rec_atoms:
        key = (a["chain_id"], a["res_seq"], a["res_name"])
        res_groups.setdefault(key, []).append(np.array([a["x"], a["y"], a["z"]], dtype=float))

    pocket_residues = []
    for (chain_id, res_seq, res_name), atom_pts in sorted(res_groups.items(), key=lambda x: x[0][1]):
        r_pts = np.vstack(atom_pts)
        d_mat = np.linalg.norm(r_pts[:, None, :] - lig_coords[None, :, :], axis=-1)
        min_d = float(d_mat.min())
        if min_d <= pocket_cutoff_A:
            pocket_residues.append({
                "residue_id": f"{res_name.capitalize()}{res_seq}",
                "res_name": res_name,
                "chain_id": chain_id,
                "res_seq": res_seq,
                "chemical_class": RESIDUE_CLASS_MAP.get(res_name, "Hydrophobic"),
                "min_distance_A": round(min_d, 3),
                "heavy_atom_count": len(atom_pts),
            })

    min_heavy_pts = []
    for line in min_rec_pdb.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("ATOM  ") and len(line) >= 54:
            el = line[76:78].strip().upper() if len(line) >= 78 else line[12:16].strip()[:1]
            if el not in ("H", "D"):
                try:
                    min_heavy_pts.append([
                        float(line[30:38].strip()),
                        float(line[38:46].strip()),
                        float(line[46:54].strip()),
                    ])
                except ValueError:
                    pass
    if min_heavy_pts:
        m_arr = np.array(min_heavy_pts, dtype=float)
        post_min_closest_A = float(np.linalg.norm(m_arr[:, None, :] - lig_coords[None, :, :], axis=-1).min())
    else:
        post_min_closest_A = float(min((r["min_distance_A"] for r in pocket_residues), default=3.0))

    min_xyz = lig_coords.min(axis=0)
    max_xyz = lig_coords.max(axis=0)
    center = 0.5 * (min_xyz + max_xyz)
    span = max_xyz - min_xyz
    size = np.minimum(s_max_A, np.maximum(s_min_A, span + 2.0 * box_padding_A))
    volume = float(size[0] * size[1] * size[2])

    return {
        "pocket_residues": pocket_residues,
        "pocket_residue_count": len(pocket_residues),
        "post_minimization_closest_heavy_atom_A": round(post_min_closest_A, 3),
        "docking_box": {
            "center_x": round(float(center[0]), 3),
            "center_y": round(float(center[1]), 3),
            "center_z": round(float(center[2]), 3),
            "size_x": round(float(size[0]), 1),
            "size_y": round(float(size[1]), 1),
            "size_z": round(float(size[2]), 1),
            "padding_A": box_padding_A,
            "volume_A3": round(volume, 1),
        },
    }


def run_stage1_target_prep(
    pdb_id: str,
    out_dir: Path,
    local_pdb: Optional[Path] = None,
    preferred_chain: Optional[str] = None,
    preferred_ligand: Optional[str] = None,
    ph: float = 7.4,
    pocket_cutoff_A: float = 5.0,
    box_padding_A: float = 4.0,
    search_dirs: Optional[List[Path]] = None,
) -> dict:
    """Executes Stage 1 end-to-end for a single target and writes all manifests & PDB/CSV files."""
    out_dir = out_dir.resolve()
    s1_dir = out_dir / "01_Stage1_Target_and_Pocket"
    s1_dir.mkdir(parents=True, exist_ok=True)

    raw_pdb = fetch_or_locate_raw_pdb(pdb_id, local_pdb, out_dir, search_dirs=search_dirs)
    sel_chain, lig_res_name, chain_lines, raw_rec_atoms, lig_atoms = isolate_primary_chain_and_ligand(
        raw_pdb, preferred_chain=preferred_chain, preferred_ligand=preferred_ligand
    )

    raw_chain_pdb = s1_dir / "receptor_isolated_chain.pdb"
    raw_chain_pdb.write_text("".join(chain_lines) + "END\n", encoding="utf-8")

    cocrystal_pdb = s1_dir / "cocrystal_reference_ligand.pdb"
    with open(cocrystal_pdb, "w", encoding="utf-8") as fl:
        for idx, a in enumerate(lig_atoms, start=1):
            fl.write(
                format_strict_pdb_hetatm(
                    idx, a["atom_name"], a["x"], a["y"], a["z"],
                    a["element"], res_name=lig_res_name, chain_id="L", res_seq=1
                )
            )
        fl.write("END\n")

    s2_dir = out_dir / "02_Stage2_Redocking_Gate"
    s2_dir.mkdir(parents=True, exist_ok=True)
    (s2_dir / "cocrystal_reference_ligand.pdb").write_bytes(cocrystal_pdb.read_bytes())

    minimized_pdb = s1_dir / "receptor_prepared_pH74_minimized.pdb"
    min_info = repair_and_minimize_receptor(raw_chain_pdb, minimized_pdb, ph=ph, max_iters=200)

    receptor_pdbqt = s1_dir / "receptor_prepared.pdbqt"
    pdbqt_compiler = prepare_receptor_pdbqt(minimized_pdb, receptor_pdbqt)

    pocket_data = compute_pocket_and_box(
        raw_rec_atoms,
        minimized_pdb,
        lig_atoms,
        pocket_cutoff_A=pocket_cutoff_A,
        box_padding_A=box_padding_A,
    )

    pocket_res_nums = {r["res_seq"] for r in pocket_data["pocket_residues"]}
    shell_lines = []
    for line in minimized_pdb.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("ATOM  ") and len(line) >= 26:
            try:
                rseq = int(line[22:26].strip())
                if rseq in pocket_res_nums:
                    shell_lines.append(line + "\n")
            except ValueError:
                pass
    shell_pdb = s1_dir / "binding_pocket_5A_shell.pdb"
    shell_pdb.write_text("".join(shell_lines) + "END\n", encoding="utf-8")

    pocket_csv = s1_dir / "pocket_residues_table.csv"
    with open(pocket_csv, "w", newline="", encoding="utf-8") as fc:
        writer = csv.DictWriter(
            fc,
            fieldnames=["residue_id", "res_name", "chain_id", "res_seq", "chemical_class", "min_distance_A", "heavy_atom_count"],
        )
        writer.writeheader()
        writer.writerows(pocket_data["pocket_residues"])

    box = pocket_data["docking_box"]
    vina_cfg = s1_dir / "vina_grid_box_config.txt"
    vina_cfg.write_text(
        f"receptor = {receptor_pdbqt.name}\n"
        f"center_x = {box['center_x']:.3f}\n"
        f"center_y = {box['center_y']:.3f}\n"
        f"center_z = {box['center_z']:.3f}\n"
        f"size_x = {box['size_x']:.1f}\n"
        f"size_y = {box['size_y']:.1f}\n"
        f"size_z = {box['size_z']:.1f}\n"
        f"exhaustiveness = 8\n",
        encoding="utf-8",
    )

    manifest = {
        "stage": 1,
        "status": "COMPLETED",
        "pdb_id": pdb_id.upper(),
        "selected_chain": sel_chain,
        "cocrystal_ligand_id": lig_res_name,
        "ligand_heavy_atoms": len(lig_atoms),
        "receptor_heavy_atoms": len(raw_rec_atoms),
        "pocket_cutoff_A": pocket_cutoff_A,
        "pocket_residue_count": pocket_data["pocket_residue_count"],
        "pocket_residues": pocket_data["pocket_residues"],
        "docking_box": box,
        "post_minimization_closest_heavy_atom_A": pocket_data["post_minimization_closest_heavy_atom_A"],
        "minimization_report": min_info,
        "pdbqt_compiler": pdbqt_compiler,
        "files": {
            "raw_pdb": str(raw_pdb),
            "minimized_receptor_pdb": str(minimized_pdb),
            "receptor_pdbqt": str(receptor_pdbqt),
            "cocrystal_ligand_pdb": str(cocrystal_pdb),
            "pocket_shell_pdb": str(shell_pdb),
            "pocket_residues_csv": str(pocket_csv),
            "vina_config_txt": str(vina_cfg),
            "minimized_receptor_sha256": sha256_file(minimized_pdb),
        },
    }

    manifest_path = out_dir / "stage1_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="AnCore v1.2 Stage 1: Multi-Target Receptor Curation & 5.0 A Pocket Mapping")
    parser.add_argument("--pdb-id", type=str, default=None, help="Single 4-character PDB ID (omit to run all 5 benchmark targets)")
    parser.add_argument("--run-all", action="store_true", help="Run all 5 benchmark targets (1IEP, 3ERT, 1UYG, 1M17, 1STP)")
    parser.add_argument("--include-4wkq", action="store_true", help="Also run the 6th custom target 4WKQ alongside the 5 benchmarks")
    parser.add_argument("--local-pdb", type=str, default=None, help="Optional path to a local .pdb file (for single-target mode)")
    parser.add_argument("--out-dir", type=str, default=None, help="Output base directory (default: workspaces/)")
    parser.add_argument("--chain", type=str, default=None, help="Preferred receptor chain ID")
    parser.add_argument("--ligand", type=str, default="AUTO", help="Preferred co-crystal ligand 3-letter code")
    parser.add_argument("--ph", type=float, default=7.4, help="Protonation pH (default: 7.4)")
    parser.add_argument("--padding", type=float, default=4.0, help="Search box padding in Angstroms (default: 4.0)")
    args, _ = parser.parse_known_args()

    cwd = Path.cwd()
    search_dirs = [
        cwd / "benchmarks",
        cwd.parent / "benchmarks",
        cwd.parent / "AnCore-Scientific-Studio" / "benchmarks",
    ]

    # Default to running all 5 benchmark targets if --pdb-id is not specified
    if args.run_all or args.pdb_id is None:
        targets = list(BENCHMARK_SUITE)
        if args.include_4wkq:
            targets.append(CUSTOM_4WKQ_TARGET)
        base_out = Path(args.out_dir) if args.out_dir else Path("workspaces")

        print("=" * 105)
        print("ANCORE v1.2 — STAGE 1 MULTI-TARGET SUITE: PDBFIXER + OPENMM AMBER14 + 5.0 Å POCKET MAPPING")
        print("=" * 105)
        print(f"{'PDB':<6} {'Ligand':<7} {'Atoms':<7} {'5Å Res':<8} {'Box Size (Å)':<18} {'Box Vol (Å³)':<14} {'Min Dist':<10} {'E_post (kcal/mol)':<18} {'Gate'}")
        print("-" * 105)

        summary_rows = []
        for t in targets:
            pid = t["pdb_id"]
            t_dir = base_out / pid
            m = run_stage1_target_prep(
                pdb_id=pid,
                out_dir=t_dir,
                preferred_chain=t["preferred_chain"],
                preferred_ligand=t["preferred_ligand"],
                ph=args.ph,
                box_padding_A=args.padding,
                search_dirs=search_dirs,
            )
            box = m["docking_box"]
            min_rep = m["minimization_report"]
            bsize = f"{box['size_x']:.1f}x{box['size_y']:.1f}x{box['size_z']:.1f}"
            min_d = m["post_minimization_closest_heavy_atom_A"]
            e_post = min_rep.get("energy_post_kcal_mol", 0.0)
            gate = "PASS" if min_d >= 2.10 else "CHECK"
            print(
                f"{pid:<6} {m['cocrystal_ligand_id']:<7} {m['ligand_heavy_atoms']:<7d} "
                f"{m['pocket_residue_count']:<8d} {bsize:<18} {box['volume_A3']:<14,.1f} "
                f"{min_d:<10.3f} {e_post:<18,.2f} {gate}"
            )
            summary_rows.append(m)

        print("=" * 105)
        suite_json = base_out / "stage1_suite_summary.json"
        suite_json.write_text(json.dumps(summary_rows, indent=2), encoding="utf-8")
        print(f"Saved multi-target Stage 1 summary -> {suite_json}")
        return

    # Single target execution
    out_dir = Path(args.out_dir) if args.out_dir else Path("workspaces") / args.pdb_id.upper()
    local_pdb = Path(args.local_pdb) if args.local_pdb else None
    m = run_stage1_target_prep(
        pdb_id=args.pdb_id,
        out_dir=out_dir,
        local_pdb=local_pdb,
        preferred_chain=args.chain,
        preferred_ligand=args.ligand,
        ph=args.ph,
        box_padding_A=args.padding,
        search_dirs=search_dirs,
    )
    box = m["docking_box"]
    min_rep = m["minimization_report"]
    hotspots_str = ", ".join(r["residue_id"] for r in m["pocket_residues"][:8])
    print("=" * 86)
    print(f"ANCORE v1.2 — STAGE 1 COMPLETED: {m['pdb_id']} (Chain {m['selected_chain']} | Ligand {m['cocrystal_ligand_id']})")
    print("=" * 86)
    print(f"  Receptor Heavy Atoms     : {m['receptor_heavy_atoms']}")
    print(f"  Ligand Heavy Atoms       : {m['ligand_heavy_atoms']} ({m['cocrystal_ligand_id']})")
    print(f"  PDBFixer + OpenMM AMBER14: E_pre = {min_rep.get('energy_pre_kcal_mol')} -> E_post = {min_rep.get('energy_post_kcal_mol')} kcal/mol")
    print(f"  Closest Heavy-Atom Dist  : {m['post_minimization_closest_heavy_atom_A']:.3f} Å (Gate >= 2.10 Å: {'PASS' if m['post_minimization_closest_heavy_atom_A'] >= 2.10 else 'CHECK'})")
    print(f"  5.0 Å Pocket Residues    : {m['pocket_residue_count']} residues ({hotspots_str}, ...)")
    print(f"  Docking Box Center (Å)   : ({box['center_x']:.2f}, {box['center_y']:.2f}, {box['center_z']:.2f})")
    print(f"  Docking Box Size (Å)     : {box['size_x']:.1f} x {box['size_y']:.1f} x {box['size_z']:.1f} (Volume = {box['volume_A3']:,.1f} Å³)")
    print(f"  Manifest Written         : {out_dir / 'stage1_manifest.json'}")
    print("=" * 86)


if __name__ == "__main__":
    main()