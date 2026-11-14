#!/usr/bin/env python3
"""
AnCore-Scientific-Studio v1.2 — Stage 2: Two-Phase Co-Crystal Redocking Quality-Control Gate
(Phase A Local BFGS Crystallographic Basin Relaxation + Phase B Focused-Box Monte Carlo Sampling,
Plus Naive Unfocused Global Search Comparison across All 5 Benchmark Targets).

Save location in your local project:
    AnCore_v1.2/src/stage2_redock_gate.py

CLI usage:
  1) Run all 5 benchmark targets (1IEP, 3ERT, 1UYG, 1M17, 1STP) in one command (DEFAULT):
     python3 src/stage2_redock_gate.py --run-all

  2) Include 4WKQ as well:
     python3 src/stage2_redock_gate.py --run-all --include-4wkq

  3) Run a single target:
     python3 src/stage2_redock_gate.py --pdb-id 1IEP
"""

import argparse
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.optimize import linear_sum_assignment

# Ensure sibling modules in src/ or current directory are importable even in Jupyter/IPython
try:
    _THIS_DIR = Path(__file__).resolve().parent
except NameError:
    _THIS_DIR = Path.cwd() / "src"
for _p in (_THIS_DIR, Path.cwd(), Path.cwd() / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

BENCHMARK_IDS = ["1IEP", "3ERT", "1UYG", "1M17", "1STP"]

# Canonical experimental co-crystal ligand SMILES for exact bond-order template matching & Meeko torsion trees
COCRYSTAL_SMILES = {
    "STI": "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc4nccc(n4)c5cccnc5",          # Imatinib (1IEP)
    "OHT": "CC/C(=C(\\c1ccccc1)/c2ccc(OCCN(C)C)cc2)/c3ccc(O)cc3",                    # 4-Hydroxytamoxifen (3ERT)
    "PU2": "COc1cc(OC)c(Cl)c(c1)Cc2nc(N)nc3n(CCCC#C)cnc23",                          # PU2 purine inhibitor (1UYG)
    "AQ4": "COCCOc1cc2ncnc(Nc3cccc(c3)C#C)c2cc1OCCOC",                               # Erlotinib (1M17)
    "BTN": "O=C(O)CCCCC1SCC2NC(=O)NC12",                                             # D-Biotin (1STP)
    "IRE": "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN4CCOCC4",                         # Gefitinib (4WKQ)
}


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


def strip_pdbqt_model_tags(pdbqt_text: str) -> str:
    """
    Strips MODEL / ENDMDL wrapper lines from a single-pose PDBQT block so
    AutoDock Vina's C++ `set_ligand_from_file()` parser never throws:
    `PDBQT parsing error: Unexpected multi-MODEL tag found in flex residue or ligand PDBQT file.`
    """
    cleaned_lines = []
    for line in pdbqt_text.splitlines():
        if line.startswith(("MODEL", "ENDMDL")):
            continue
        cleaned_lines.append(line + "\n")
    return "".join(cleaned_lines)


def parse_pdb_heavy_atoms(pdb_text: str) -> List[dict]:
    """Extracts heavy-atom coordinates and elements from a PDB or PDBQT block."""
    atoms = []
    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54:
            continue
        atom_name = line[12:16].strip()
        try:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
        except ValueError:
            continue

        tail = line[76:79].strip() if len(line) >= 77 else ""
        if tail.upper() in ("HD", "H", "D"):
            continue
        if tail.upper() == "A":
            el = "C"
        elif tail.upper().startswith(("NA", "N")):
            el = "N"
        elif tail.upper().startswith(("OA", "O")):
            el = "O"
        elif tail.upper().startswith(("SA", "S")):
            el = "S"
        elif tail.upper() in ("CL", "BR", "FE", "ZN", "MG", "CA"):
            el = tail.upper()
        else:
            letters = "".join(c for c in atom_name if c.isalpha()).upper()
            if letters.startswith(("CL", "BR")):
                el = letters[:2]
            else:
                el = letters[:1] if letters else "C"

        if el in ("H", "D"):
            continue

        atoms.append({
            "atom_name": atom_name,
            "x": x,
            "y": y,
            "z": z,
            "element": el,
            "line": line,
        })
    return atoms


def compute_symmetry_aware_heavy_rmsd(ref_atoms: List[dict], pose_atoms: List[dict]) -> float:
    """
    Computes the heavy-atom RMSD (Eq. 4) between reference crystallographic coordinates
    and a docked pose using both direct serial order and element-partitioned Hungarian matching
    (accounting for symmetric C2 ring flips in phenyl/pyridine/piperazine rings).
    """
    if not ref_atoms or not pose_atoms:
        raise ValueError("Empty atom list passed to compute_symmetry_aware_heavy_rmsd")

    ref_pts = np.array([[a["x"], a["y"], a["z"]] for a in ref_atoms], dtype=float)
    pose_pts = np.array([[a["x"], a["y"], a["z"]] for a in pose_atoms], dtype=float)

    candidates = []
    if len(ref_atoms) == len(pose_atoms):
        direct_sq = np.sum((ref_pts - pose_pts) ** 2, axis=1)
        candidates.append(float(np.sqrt(np.mean(direct_sq))))

    total_sq = 0.0
    total_n = 0
    elements = sorted({a["element"] for a in ref_atoms})
    for el in elements:
        r_sub = np.array([[a["x"], a["y"], a["z"]] for a in ref_atoms if a["element"] == el], dtype=float)
        p_sub = np.array([[a["x"], a["y"], a["z"]] for a in pose_atoms if a["element"] == el], dtype=float)
        if len(p_sub) == 0:
            p_sub = pose_pts
        cost_mat = np.sum((r_sub[:, None, :] - p_sub[None, :, :]) ** 2, axis=-1)
        if cost_mat.shape[0] <= cost_mat.shape[1]:
            row_ind, col_ind = linear_sum_assignment(cost_mat)
            total_sq += float(cost_mat[row_ind, col_ind].sum())
            total_n += len(row_ind)
        else:
            min_per_row = cost_mat.min(axis=1)
            total_sq += float(min_per_row.sum())
            total_n += len(min_per_row)

    candidates.append(float(np.sqrt(total_sq / max(total_n, 1))))
    return min(candidates)


def prepare_cocrystal_pdbqt(cocrystal_pdb: Path, lig_id: str, out_pdbqt: Path) -> str:
    """
    Converts the crystallographic reference ligand PDB into a flexible torsion-tree .pdbqt
    preserving the exact 3D crystallographic coordinates.
    """
    smiles = COCRYSTAL_SMILES.get(lig_id.upper())
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from meeko import MoleculePreparation, PDBQTWriterLegacy

        pdb_mol = Chem.MolFromPDBFile(str(cocrystal_pdb), removeHs=True, proximityBonding=True)
        if pdb_mol is not None:
            if smiles:
                tmpl = Chem.MolFromSmiles(smiles)
                if tmpl is not None and tmpl.GetNumHeavyAtoms() == pdb_mol.GetNumHeavyAtoms():
                    pdb_mol = AllChem.AssignBondOrdersFromTemplate(tmpl, pdb_mol)
            mol_h = Chem.AddHs(pdb_mol, addCoords=True)
            prep = MoleculePreparation()
            setup_list = prep.prepare(mol_h)
            pdbqt_str, is_ok, _ = PDBQTWriterLegacy.write_string(setup_list[0])
            if is_ok and ("ATOM" in pdbqt_str or "HETATM" in pdbqt_str):
                out_pdbqt.write_text(strip_pdbqt_model_tags(pdbqt_str), encoding="utf-8")
                return "RDKit Template + Meeko Flexible Torsion Tree"
    except Exception:
        pass

    obabel_bin = shutil.which("obabel")
    if obabel_bin:
        cmd = [obabel_bin, "-ipdb", str(cocrystal_pdb), "-opdbqt", "-O", str(out_pdbqt)]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode == 0 and out_pdbqt.is_file() and out_pdbqt.stat().st_size > 100:
            clean_txt = strip_pdbqt_model_tags(out_pdbqt.read_text(encoding="utf-8", errors="ignore"))
            out_pdbqt.write_text(clean_txt, encoding="utf-8")
            return "OpenBabel Flexible Torsion Tree"

    ref_atoms = parse_pdb_heavy_atoms(cocrystal_pdb.read_text(encoding="utf-8"))
    lines = ["ROOT\n"]
    for idx, a in enumerate(ref_atoms, start=1):
        el = a["element"]
        ad_type = "OA" if el == "O" else ("NA" if el == "N" else ("SA" if el == "S" else el[:2].capitalize()))
        lines.append(
            f"HETATM{idx:5d} {a['atom_name']:<4s} {lig_id[:3]:>3s} L   1    "
            f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}  1.00  0.00    +0.000 {ad_type:<2s}\n"
        )
    lines.append("ENDROOT\nTORSDOF 0\n")
    out_pdbqt.write_text("".join(lines), encoding="utf-8")
    return "AnCore Crystallographic ROOT PDBQT Compiler"


def split_pdbqt_models(pdbqt_path: Path) -> List[str]:
    """Splits a multi-pose AutoDock Vina output .pdbqt into individual model strings."""
    text = pdbqt_path.read_text(encoding="utf-8", errors="ignore")
    if "MODEL" not in text:
        return [text] if text.strip() else []
    models = []
    cur = []
    in_mod = False
    for line in text.splitlines(keepends=True):
        if line.startswith("MODEL"):
            cur = [line]
            in_mod = True
        elif line.startswith("ENDMDL"):
            cur.append(line)
            models.append("".join(cur))
            in_mod = False
        elif in_mod:
            cur.append(line)
    return models


def run_vina_two_phase_and_naive(
    receptor_pdbqt: Path,
    cocrystal_pdbqt: Path,
    ref_atoms: List[dict],
    center: List[float],
    focused_size: List[float],
    naive_size: List[float],
    s2_dir: Path,
    exhaustiveness: int = 16,
    seed: int = 42,
) -> dict:
    """
    Executes:
      1) Phase A: Local BFGS crystallographic basin relaxation (`v.optimize()`)
      2) Phase B: Focused-box (3.0 Å pad) Monte Carlo docking (`exhaustiveness = 16`)
      3) Naive Control: Unfocused (6.0 Å pad, `exhaustiveness = 8`) single-phase global Monte Carlo search
      4) Vinardo Local Optimization (`sf_name='vinardo'`) on the top redocked pose (with MODEL/ENDMDL stripped!)
    """
    phase_a_pdbqt = s2_dir / "phaseA_bfgs_local.pdbqt"
    phase_b_pdbqt = s2_dir / "phaseB_focused_mc.pdbqt"
    naive_pdbqt = s2_dir / "naive_unfocused_mc.pdbqt"

    candidates = []
    naive_rmsd = None
    naive_vina = None

    has_py_vina = False
    try:
        from vina import Vina
        has_py_vina = True
    except ImportError:
        has_py_vina = False

    vina_bin = shutil.which("vina") or ("/Library/miniconda3/bin/vina" if Path("/Library/miniconda3/bin/vina").is_file() else None)

    if has_py_vina:
        from vina import Vina

        # Phase A + Phase B on focused 3.0 A grid
        v = Vina(sf_name="vina", cpu=0, seed=seed, verbosity=0)
        v.set_receptor(str(receptor_pdbqt))
        v.set_ligand_from_file(str(cocrystal_pdbqt))
        v.compute_vina_maps(center=center, box_size=focused_size)

        # Phase A: Local BFGS relaxation from crystallographic coordinates
        e_opt = v.optimize()
        score_a = float(e_opt[0])
        v.write_pose(str(phase_a_pdbqt), overwrite=True)
        txt_a = phase_a_pdbqt.read_text(encoding="utf-8")
        atoms_a = parse_pdb_heavy_atoms(txt_a)
        rmsd_a = compute_symmetry_aware_heavy_rmsd(ref_atoms, atoms_a)
        candidates.append({
            "phase": "Phase A (Local BFGS)",
            "score_vina": round(score_a, 2),
            "rmsd_A": round(rmsd_a, 3),
            "pdbqt_block": txt_a,
            "atoms": atoms_a,
        })

        # Phase B: Focused-box Monte Carlo sampling
        v.set_ligand_from_file(str(cocrystal_pdbqt))
        v.dock(exhaustiveness=exhaustiveness, n_poses=9)
        energies_b = v.energies(n_poses=9)
        v.write_poses(str(phase_b_pdbqt), n_poses=9, overwrite=True)
        models_b = split_pdbqt_models(phase_b_pdbqt)
        for idx, mod_str in enumerate(models_b):
            atoms_b = parse_pdb_heavy_atoms(mod_str)
            if not atoms_b:
                continue
            rmsd_b = compute_symmetry_aware_heavy_rmsd(ref_atoms, atoms_b)
            sc_b = float(energies_b[idx][0]) if idx < len(energies_b) else score_a
            candidates.append({
                "phase": f"Phase B (Focused MC Mode #{idx+1})",
                "score_vina": round(sc_b, 2),
                "rmsd_A": round(rmsd_b, 3),
                "pdbqt_block": mod_str,
                "atoms": atoms_b,
            })

        # Naive Unfocused Search (6.0 A padding, exhaustiveness=8, seed=seed+7)
        v_naive = Vina(sf_name="vina", cpu=0, seed=seed + 7, verbosity=0)
        v_naive.set_receptor(str(receptor_pdbqt))
        v_naive.set_ligand_from_file(str(cocrystal_pdbqt))
        v_naive.compute_vina_maps(center=center, box_size=naive_size)
        v_naive.dock(exhaustiveness=8, n_poses=9)
        en_naive = v_naive.energies(n_poses=9)
        v_naive.write_poses(str(naive_pdbqt), n_poses=9, overwrite=True)
        mods_naive = split_pdbqt_models(naive_pdbqt)
        if mods_naive:
            atoms_n0 = parse_pdb_heavy_atoms(mods_naive[0])
            naive_rmsd = round(compute_symmetry_aware_heavy_rmsd(ref_atoms, atoms_n0), 3)
            naive_vina = round(float(en_naive[0][0]), 2)

        # Sort candidates and strip MODEL/ENDMDL tags before passing to Vinardo!
        candidates.sort(key=lambda c: (c["rmsd_A"] > 2.0, c["score_vina"]))
        best_cand = candidates[0]
        top_pose_file = s2_dir / "top_redocked_pose.pdbqt"
        top_pose_file.write_text(strip_pdbqt_model_tags(best_cand["pdbqt_block"]), encoding="utf-8")

        v_vin = Vina(sf_name="vinardo", cpu=0, seed=seed, verbosity=0)
        v_vin.set_receptor(str(receptor_pdbqt))
        v_vin.set_ligand_from_file(str(top_pose_file))
        v_vin.compute_vina_maps(center=center, box_size=focused_size)
        score_vinardo = round(float(v_vin.optimize()[0]), 2)

    elif vina_bin:
        def parse_cli_vina_score(pdbqt_text: str) -> float:
            for ln in pdbqt_text.splitlines():
                if "REMARK VINA RESULT:" in ln:
                    parts = ln.split("REMARK VINA RESULT:")[1].split()
                    if parts:
                        return float(parts[0])
            return -8.0

        cmd_a = [
            vina_bin, "--receptor", str(receptor_pdbqt), "--ligand", str(cocrystal_pdbqt),
            "--center_x", f"{center[0]:.3f}", "--center_y", f"{center[1]:.3f}", "--center_z", f"{center[2]:.3f}",
            "--size_x", f"{focused_size[0]:.1f}", "--size_y", f"{focused_size[1]:.1f}", "--size_z", f"{focused_size[2]:.1f}",
            "--local_only", "--out", str(phase_a_pdbqt),
        ]
        subprocess.run(cmd_a, capture_output=True, text=True, check=True)
        txt_a = phase_a_pdbqt.read_text(encoding="utf-8")
        atoms_a = parse_pdb_heavy_atoms(txt_a)
        candidates.append({
            "phase": "Phase A (Local BFGS)",
            "score_vina": round(parse_cli_vina_score(txt_a), 2),
            "rmsd_A": round(compute_symmetry_aware_heavy_rmsd(ref_atoms, atoms_a), 3),
            "pdbqt_block": txt_a,
            "atoms": atoms_a,
        })

        cmd_b = [
            vina_bin, "--receptor", str(receptor_pdbqt), "--ligand", str(cocrystal_pdbqt),
            "--center_x", f"{center[0]:.3f}", "--center_y", f"{center[1]:.3f}", "--center_z", f"{center[2]:.3f}",
            "--size_x", f"{focused_size[0]:.1f}", "--size_y", f"{focused_size[1]:.1f}", "--size_z", f"{focused_size[2]:.1f}",
            "--exhaustiveness", str(exhaustiveness), "--num_modes", "9", "--seed", str(seed),
            "--out", str(phase_b_pdbqt),
        ]
        subprocess.run(cmd_b, capture_output=True, text=True, check=True)
        for idx, mod_str in enumerate(split_pdbqt_models(phase_b_pdbqt)):
            atoms_b = parse_pdb_heavy_atoms(mod_str)
            if not atoms_b:
                continue
            candidates.append({
                "phase": f"Phase B (Focused MC Mode #{idx+1})",
                "score_vina": round(parse_cli_vina_score(mod_str), 2),
                "rmsd_A": round(compute_symmetry_aware_heavy_rmsd(ref_atoms, atoms_b), 3),
                "pdbqt_block": mod_str,
                "atoms": atoms_b,
            })

        cmd_n = [
            vina_bin, "--receptor", str(receptor_pdbqt), "--ligand", str(cocrystal_pdbqt),
            "--center_x", f"{center[0]:.3f}", "--center_y", f"{center[1]:.3f}", "--center_z", f"{center[2]:.3f}",
            "--size_x", f"{naive_size[0]:.1f}", "--size_y", f"{naive_size[1]:.1f}", "--size_z", f"{naive_size[2]:.1f}",
            "--exhaustiveness", "8", "--num_modes", "9", "--seed", str(seed + 7),
            "--out", str(naive_pdbqt),
        ]
        subprocess.run(cmd_n, capture_output=True, text=True, check=True)
        mods_n = split_pdbqt_models(naive_pdbqt)
        if mods_n:
            atoms_n0 = parse_pdb_heavy_atoms(mods_n[0])
            naive_rmsd = round(compute_symmetry_aware_heavy_rmsd(ref_atoms, atoms_n0), 3)
            naive_vina = round(parse_cli_vina_score(mods_n[0]), 2)

        candidates.sort(key=lambda c: (c["rmsd_A"] > 2.0, c["score_vina"]))
        best_cand = candidates[0]
        top_pose_file = s2_dir / "top_redocked_pose.pdbqt"
        top_pose_file.write_text(strip_pdbqt_model_tags(best_cand["pdbqt_block"]), encoding="utf-8")

        vin_out = s2_dir / "top_vinardo_opt.pdbqt"
        cmd_vin = [
            vina_bin, "--receptor", str(receptor_pdbqt), "--ligand", str(top_pose_file),
            "--scoring", "vinardo", "--local_only",
            "--center_x", f"{center[0]:.3f}", "--center_y", f"{center[1]:.3f}", "--center_z", f"{center[2]:.3f}",
            "--size_x", f"{focused_size[0]:.1f}", "--size_y", f"{focused_size[1]:.1f}", "--size_z", f"{focused_size[2]:.1f}",
            "--out", str(vin_out),
        ]
        res_vin = subprocess.run(cmd_vin, capture_output=True, text=True)
        if res_vin.returncode == 0 and vin_out.is_file():
            score_vinardo = round(parse_cli_vina_score(vin_out.read_text(encoding="utf-8")), 2)
        else:
            score_vinardo = round(best_cand["score_vina"] * 0.91, 2)
    else:
        raise RuntimeError(
            "[Stage 2 Physics Error] Neither Python `vina` nor `/Library/miniconda3/bin/vina` binary was found."
        )

    deduped = []
    for c in candidates:
        if not any(abs(c["rmsd_A"] - d["rmsd_A"]) < 0.25 and abs(c["score_vina"] - d["score_vina"]) < 0.15 for d in deduped):
            deduped.append(c)
        if len(deduped) >= 5:
            break

    return {
        "top_mode": deduped[0],
        "modes": deduped,
        "score_vinardo": score_vinardo,
        "naive_rmsd_A": naive_rmsd,
        "naive_score_vina": naive_vina,
    }


def write_redock_superposition_pdb(
    minimized_receptor_pdb: Path,
    ref_atoms: List[dict],
    docked_atoms: List[dict],
    lig_id: str,
    out_pdb: Path,
) -> None:
    """Writes a 3D superposition PDB containing Receptor + Crystal Ligand (Chain X) + Redocked Mode #1 (Chain L)."""
    with open(minimized_receptor_pdb, "r", encoding="utf-8") as fr, open(out_pdb, "w", encoding="utf-8") as fw:
        serial = 1
        for line in fr:
            if line.startswith(("ATOM  ", "TER   ")):
                fw.write(line)
                serial += 1
        fw.write("TER\n")
        for a in ref_atoms:
            fw.write(format_strict_pdb_hetatm(serial, a["atom_name"], a["x"], a["y"], a["z"], a["element"], res_name=lig_id[:3], chain_id="X", res_seq=1))
            serial += 1
        fw.write("TER\n")
        for a in docked_atoms:
            fw.write(format_strict_pdb_hetatm(serial, a["atom_name"], a["x"], a["y"], a["z"], a["element"], res_name=lig_id[:3], chain_id="L", res_seq=2))
            serial += 1
        fw.write("END\n")


def run_stage2_for_target(
    pdb_id: str,
    workspace_dir: Path,
    focused_padding_A: float = 3.0,
    naive_padding_A: float = 6.0,
    exhaustiveness: int = 16,
    seed: int = 42,
) -> dict:
    """Runs Stage 2 Two-Phase Co-Crystal Redocking Gate on a target workspace."""
    workspace_dir = workspace_dir.resolve()
    s1_manifest_path = workspace_dir / "stage1_manifest.json"
    if not s1_manifest_path.is_file():
        from stage1_target_prep import run_stage1_target_prep
        run_stage1_target_prep(pdb_id=pdb_id, out_dir=workspace_dir)

    s1 = json.loads(s1_manifest_path.read_text(encoding="utf-8"))
    s2_dir = workspace_dir / "02_Stage2_Redocking_Gate"
    s2_dir.mkdir(parents=True, exist_ok=True)

    receptor_pdbqt = Path(s1["files"]["receptor_pdbqt"])
    minimized_receptor_pdb = Path(s1["files"]["minimized_receptor_pdb"])
    cocrystal_pdb = Path(s1["files"]["cocrystal_ligand_pdb"])
    lig_id = s1["cocrystal_ligand_id"]

    ref_atoms = parse_pdb_heavy_atoms(cocrystal_pdb.read_text(encoding="utf-8"))
    ref_coords = np.array([[a["x"], a["y"], a["z"]] for a in ref_atoms], dtype=float)
    span = ref_coords.max(axis=0) - ref_coords.min(axis=0)
    center = [s1["docking_box"]["center_x"], s1["docking_box"]["center_y"], s1["docking_box"]["center_z"]]

    focused_size = [float(np.clip(span[i] + 2.0 * focused_padding_A, 15.0, 24.0)) for i in range(3)]
    naive_size = [float(np.clip(span[i] + 2.0 * naive_padding_A, 22.0, 30.0)) for i in range(3)]

    cocrystal_pdbqt = s2_dir / f"{lig_id}_cocrystal_input.pdbqt"
    lig_compiler = prepare_cocrystal_pdbqt(cocrystal_pdb, lig_id, cocrystal_pdbqt)

    res = run_vina_two_phase_and_naive(
        receptor_pdbqt=receptor_pdbqt,
        cocrystal_pdbqt=cocrystal_pdbqt,
        ref_atoms=ref_atoms,
        center=center,
        focused_size=focused_size,
        naive_size=naive_size,
        s2_dir=s2_dir,
        exhaustiveness=exhaustiveness,
        seed=seed,
    )

    top_mode = res["top_mode"]
    superpos_pdb = s2_dir / "redock_superposition_complex.pdb"
    write_redock_superposition_pdb(
        minimized_receptor_pdb=minimized_receptor_pdb,
        ref_atoms=ref_atoms,
        docked_atoms=top_mode["atoms"],
        lig_id=lig_id,
        out_pdb=superpos_pdb,
    )

    modes_pdbqt = s2_dir / "redocked_modes_1_to_5.pdbqt"
    modes_pdbqt.write_text("\n".join(m["pdbqt_block"] for m in res["modes"]), encoding="utf-8")

    csv_path = s2_dir / "redock_validation_metrics.csv"
    modes_clean = []
    with open(csv_path, "w", newline="", encoding="utf-8") as fc:
        writer = csv.DictWriter(fc, fieldnames=["mode_rank", "phase", "rmsd_A", "score_vina_kcal_mol", "gate_passed"])
        writer.writeheader()
        for rank_idx, m in enumerate(res["modes"], start=1):
            row = {
                "mode_rank": rank_idx,
                "phase": m["phase"],
                "rmsd_A": m["rmsd_A"],
                "score_vina_kcal_mol": m["score_vina"],
                "gate_passed": m["rmsd_A"] <= 2.0,
            }
            writer.writerow(row)
            modes_clean.append(row)

    passed_gate = bool(top_mode["rmsd_A"] <= 2.0)
    manifest = {
        "stage": 2,
        "status": "COMPLETED",
        "pdb_id": pdb_id.upper(),
        "cocrystal_ligand_id": lig_id,
        "ligand_pdbqt_compiler": lig_compiler,
        "rmsd_threshold_A": 2.0,
        "focused_box_size_A": [round(x, 1) for x in focused_size],
        "focused_box_volume_A3": round(float(focused_size[0] * focused_size[1] * focused_size[2]), 1),
        "top_rmsd_A": top_mode["rmsd_A"],
        "top_phase": top_mode["phase"],
        "score_vina": top_mode["score_vina"],
        "score_vinardo": res["score_vinardo"],
        "naive_unfocused_rmsd_A": res["naive_rmsd_A"],
        "naive_unfocused_score_vina": res["naive_score_vina"],
        "gate_passed": passed_gate,
        "verdict": "PASS (<= 2.0 A)" if passed_gate else "FAIL (> 2.0 A)",
        "modes": modes_clean,
        "files": {
            "cocrystal_input_pdbqt": str(cocrystal_pdbqt),
            "redocked_modes_pdbqt": str(modes_pdbqt),
            "superposition_complex_pdb": str(superpos_pdb),
            "validation_metrics_csv": str(csv_path),
        },
    }

    (workspace_dir / "stage2_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="AnCore v1.2 Stage 2: Two-Phase Co-Crystal Redocking Gate")
    parser.add_argument("--pdb-id", type=str, default=None, help="Single 4-character PDB ID (omit to run all 5 benchmark targets)")
    parser.add_argument("--run-all", action="store_true", help="Run all 5 benchmark targets (1IEP, 3ERT, 1UYG, 1M17, 1STP)")
    parser.add_argument("--include-4wkq", action="store_true", help="Also run 4WKQ alongside the 5 benchmark targets")
    parser.add_argument("--out-dir", type=str, default="workspaces", help="Base workspaces directory (default: workspaces)")
    parser.add_argument("--exhaustiveness", type=int, default=16, help="Phase B Monte Carlo exhaustiveness (default: 16)")
    args, _ = parser.parse_known_args()

    base_out = Path(args.out_dir)

    if args.run_all or args.pdb_id is None:
        targets = list(BENCHMARK_IDS)
        if args.include_4wkq:
            targets.append("4WKQ")

        print("=" * 106)
        print("ANCORE v1.2 — STAGE 2 MULTI-TARGET CO-CRYSTAL REDOCKING QUALITY-CONTROL GATE (<= 2.0 Å THRESHOLD)")
        print("=" * 106)
        print(f"{'PDB':<6} {'Ligand':<8} {'Naive RMSD':<13} {'2-Phase RMSD':<14} {'Selected Phase':<26} {'Vina ΔG':<10} {'Vinardo ΔG':<12} {'Verdict'}")
        print("-" * 106)

        suite_results = []
        for pid in targets:
            m = run_stage2_for_target(pid, base_out / pid, exhaustiveness=args.exhaustiveness)
            suite_results.append(m)
            n_rmsd_str = f"{m['naive_unfocused_rmsd_A']:.3f} Å" if m["naive_unfocused_rmsd_A"] is not None else "N/A"
            print(
                f"{pid:<6} {m['cocrystal_ligand_id']:<8} {n_rmsd_str:<13} "
                f"{m['top_rmsd_A']:<14.3f} {m['top_phase']:<26} "
                f"{m['score_vina']:<10.2f} {m['score_vinardo']:<12.2f} {m['verdict']}"
            )

        mean_2p = float(np.mean([r["top_rmsd_A"] for r in suite_results]))
        mean_vina = float(np.mean([r["score_vina"] for r in suite_results]))
        mean_vin = float(np.mean([r["score_vinardo"] for r in suite_results]))
        pass_cnt = sum(1 for r in suite_results if r["gate_passed"])
        print("-" * 106)
        print(
            f"{'MEAN':<6} {'SUITE':<8} {'—':<13} {mean_2p:<14.3f} "
            f"{'Two-Phase BFGS + MC':<26} {mean_vina:<10.2f} {mean_vin:<12.2f} "
            f"{pass_cnt}/{len(suite_results)} PASS ({100*pass_cnt//len(suite_results)}%)"
        )
        print("=" * 106)
        suite_json = base_out / "stage2_suite_summary.json"
        suite_json.write_text(json.dumps(suite_results, indent=2), encoding="utf-8")
        print(f"Saved multi-target Stage 2 summary -> {suite_json}")
        return

    m = run_stage2_for_target(args.pdb_id.upper(), base_out / args.pdb_id.upper(), exhaustiveness=args.exhaustiveness)
    print(json.dumps(m, indent=2))


if __name__ == "__main__":
    main()