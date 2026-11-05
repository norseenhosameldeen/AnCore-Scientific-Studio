#!/usr/bin/env python3
"""
AnCore — Stage 1: Target Curation, Repair & 3D Search Box (AnCore Studio v2.2)
==============================================================================
Authors: Norseen Hossam El-Deen & Ahmed M. Senary
Affiliation: Faculty of Pharmacy, Galala University, Egypt
License: Apache-2.0 (100% Original Open-Source Implementation)

Description:
    Standalone Stage 1 module for the AnCore Native macOS Molecular Docking Workbench.
    Features an original dark-mode scientific studio interface (PySide6 + 3Dmol.js) and
    a 5-target gold-standard benchmark suite:
    1. 5-Target Gold-Standard Benchmark Suite (CASF-2016, DUD-E, LIT-PCBA):
       - 1IEP: ABL1 Tyrosine Kinase (STI / Imatinib) — DFG-out allosteric tunnel
       - 3ERT: Estrogen Receptor alpha (OHT / 4-Hydroxytamoxifen) — Steroid cavity
       - 1UYG: HSP90-alpha Chaperone (PU3 / Purine inhibitor) — Polar ATPase pocket
       - 1M17: EGFR Tyrosine Kinase (AQ4 / Erlotinib) — Open ATP hinge cleft
       - 1STP: Streptavidin (BTN / Biotin) — CASF-2016 high-affinity H-bond cage
    2. Batch & Single-Target Curation:
       - Curates all 5 benchmark targets (`pdb_id="ALL"`) or any single target.
       - RCSB PDB download, chain isolation, PDBFixer/OpenMM heavy-atom repair & pH 7.4 protonation.
       - Flexible ligand PDBQT preparation (preserving rotatable bonds for Stage 2 redocking).
       - 5.0 A binding pocket contact mapping & 3D docking box calculation.
    3. AnCore Dark Studio GUI (PySide6 + 3Dmol.js):
       - Left Benchmark Suite Navigator table + clean 2-column parameter grid.
       - Widescreen 3D molecular viewport with HUD overlay, cavity surface toggle, and toolbar.
       - Bottom tabbed table: [Benchmark Suite Overview (5 Targets)] & [5.0 A Pocket Residues].
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import os
import shutil
import subprocess
import sys
import time
import urllib.request
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("AnCore.Stage1")

# ==============================================================================
# GOLD-STANDARD 5-TARGET BENCHMARK SUITE (CASF-2016 / DUD-E / LIT-PCBA)
# ==============================================================================
BENCHMARK_TARGETS: Dict[str, Dict[str, Any]] = {
    "1IEP": {
        "name": "ABL1 Tyrosine Kinase",
        "ligand": "STI",
        "ligand_name": "Imatinib",
        "chain": "A",
        "resolution": "2.10 A",
        "benchmark": "DUD-E / Vina Tutorial",
        "pocket_type": "Deep DFG-out allosteric tunnel (Gatekeeper THR315)",
        "default_padding": 4.0,
    },
    "3ERT": {
        "name": "Estrogen Receptor alpha",
        "ligand": "OHT",
        "ligand_name": "4-Hydroxytamoxifen",
        "chain": "A",
        "resolution": "1.90 A",
        "benchmark": "DUD-E & LIT-PCBA",
        "pocket_type": "Enclosed hydrophobic steroid nuclear receptor cavity",
        "default_padding": 4.0,
    },
    "1UYG": {
        "name": "HSP90-alpha Chaperone",
        "ligand": "PU3",
        "ligand_name": "Purine Inhibitor",
        "chain": "A",
        "resolution": "1.90 A",
        "benchmark": "DUD-E Benchmark",
        "pocket_type": "Rigid polar ATP-binding Bergerat fold",
        "default_padding": 4.0,
    },
    "1M17": {
        "name": "EGFR Tyrosine Kinase",
        "ligand": "AQ4",
        "ligand_name": "Erlotinib",
        "chain": "A",
        "resolution": "2.60 A",
        "benchmark": "DUD-E Oncology",
        "pocket_type": "Open ATP hinge-binding cleft (MET769 H-bond)",
        "default_padding": 4.0,
    },
    "1STP": {
        "name": "Streptavidin",
        "ligand": "BTN",
        "ligand_name": "Biotin (Vitamin B7)",
        "chain": "A",
        "resolution": "1.40 A",
        "benchmark": "CASF-2016 Gold Standard",
        "pocket_type": "Compact femtomolar beta-barrel H-bond cage",
        "default_padding": 4.0,
    },
}

STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP", "CYX",
}

WATER_RES = {"HOH", "WAT", "TIP", "TIP3", "SOL", "DOD"}

BUFFER_IONS = {
    "SO4", "PO4", "GOL", "EDO", "DMS", "ACT", "ACE", "CL", "NA", "MG",
    "ZN", "CA", "K", "MN", "FE", "IOD", "BR", "FMT", "TRS", "PEG", "PGE",
    "MPD", "CIT", "NO3", "SCN", "CAC", "BME",
}

AD4_PROTON_MAP = {
    "C": "C", "N": "NA", "O": "OA", "S": "SA", "P": "P",
    "F": "F", "CL": "Cl", "BR": "Br", "I": "I", "H": "HD",
    "ZN": "Zn", "MG": "Mg", "CA": "Ca", "FE": "Fe", "MN": "Mn",
}


# ==============================================================================
# CORE SCIENTIFIC ENGINE
# ==============================================================================
def download_pdb(pdb_id: str, out_dir: Path) -> Path:
    """Download crystal structure from RCSB PDB or return cached file."""
    pdb_id = pdb_id.strip().upper()
    out_dir.mkdir(parents=True, exist_ok=True)
    pdb_path = out_dir / f"{pdb_id}_raw.pdb"
    if pdb_path.exists() and pdb_path.stat().st_size > 1000:
        return pdb_path

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    logger.info(f"Fetching {pdb_id} from RCSB PDB ({url})...")
    req = urllib.request.Request(url, headers={"User-Agent": "AnCore-Workbench/2.2"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        data = resp.read().decode("utf-8", errors="replace")
    if "ATOM" not in data and "HETATM" not in data:
        raise ValueError(f"Downloaded file for {pdb_id} does not contain valid ATOM records.")
    pdb_path.write_text(data, encoding="utf-8")
    return pdb_path


def isolate_receptor_and_ligand(
    raw_pdb_path: Path,
    work_dir: Path,
    target_chain: Optional[str] = "A",
    ref_ligand: Optional[str] = None,
) -> Tuple[Path, Path, str, str, List[Tuple[float, float, float]]]:
    """Separate clean protein chain from co-crystallized reference ligand."""
    lines = raw_pdb_path.read_text(encoding="utf-8", errors="replace").splitlines()

    available_chains: List[str] = []
    for ln in lines:
        if ln.startswith("ATOM  "):
            ch = ln[21:22].strip() or "A"
            if ch not in available_chains:
                available_chains.append(ch)

    chosen_chain = target_chain.strip().upper() if target_chain else (available_chains[0] if available_chains else "A")
    if chosen_chain not in available_chains and available_chains:
        chosen_chain = available_chains[0]

    receptor_lines: List[str] = []
    het_by_key: Dict[Tuple[str, str, str], List[str]] = {}

    for ln in lines:
        if ln.startswith("ATOM  "):
            ch = ln[21:22].strip() or "A"
            alt = ln[16:17].strip()
            if alt not in ("", "A", "1"):
                continue
            resname = ln[17:20].strip()
            if ch == chosen_chain and resname in STANDARD_AA:
                receptor_lines.append(ln)
        elif ln.startswith("HETATM"):
            resname = ln[17:20].strip()
            if resname in WATER_RES or resname in BUFFER_IONS:
                continue
            ch = ln[21:22].strip() or "A"
            resseq = ln[22:26].strip()
            alt = ln[16:17].strip()
            if alt not in ("", "A", "1"):
                continue
            het_by_key.setdefault((resname, ch, resseq), []).append(ln)

    if not receptor_lines:
        raise ValueError(f"No standard protein ATOM records found for chain {chosen_chain}.")

    chosen_het_key: Optional[Tuple[str, str, str]] = None
    req_lig = ref_ligand.strip().upper() if ref_ligand else None

    if req_lig:
        candidates = [k for k in het_by_key if k[0] == req_lig and k[1] == chosen_chain]
        if not candidates:
            candidates = [k for k in het_by_key if k[0] == req_lig]
        if candidates:
            chosen_het_key = max(candidates, key=lambda k: len(het_by_key[k]))

    if chosen_het_key is None and het_by_key:
        chain_cands = [k for k in het_by_key if k[1] == chosen_chain]
        pool = chain_cands if chain_cands else list(het_by_key.keys())
        chosen_het_key = max(pool, key=lambda k: len(het_by_key[k]))

    if chosen_het_key is None:
        raise ValueError(f"No co-crystallized ligand found in {raw_pdb_path.name}.")

    detected_lig = chosen_het_key[0]
    lig_lines = het_by_key[chosen_het_key]

    lig_coords: List[Tuple[float, float, float]] = []
    for ln in lig_lines:
        elem = (ln[76:78].strip() if len(ln) >= 78 else ln[12:16].strip()[:1]).upper()
        if elem == "H":
            continue
        try:
            x = float(ln[30:38])
            y = float(ln[38:46])
            z = float(ln[46:54])
            lig_coords.append((x, y, z))
        except ValueError:
            continue

    receptor_raw_pdb = work_dir / "receptor_isolated.pdb"
    receptor_raw_pdb.write_text("\n".join(receptor_lines + ["TER", "END"]) + "\n", encoding="utf-8")

    crystal_lig_pdb = work_dir / "crystal_ligand.pdb"
    crystal_lig_pdb.write_text("\n".join(lig_lines + ["END"]) + "\n", encoding="utf-8")

    return receptor_raw_pdb, crystal_lig_pdb, chosen_chain, detected_lig, lig_coords


def repair_and_protonate_receptor(isolated_pdb: Path, clean_pdb: Path, ph: float = 7.4) -> str:
    """Repair missing heavy atoms and protonate receptor at target pH using PDBFixer + OpenMM."""
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile

        fixer = PDBFixer(filename=str(isolated_pdb))
        fixer.findMissingResidues()
        fixer.missingResidues = {}
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(keepWater=False)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(ph)

        with open(clean_pdb, "w", encoding="utf-8") as fh:
            PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
        return f"PDBFixer + OpenMM (pH {ph:.1f})"
    except Exception as exc:
        logger.warning(f"PDBFixer fallback ({exc}); using isolated receptor chain.")
        shutil.copyfile(isolated_pdb, clean_pdb)
        return "Direct Chain Isolation"


def receptor_pdb_to_pdbqt(clean_pdb: Path, out_pdbqt: Path) -> str:
    """Convert clean receptor PDB to rigid AutoDock PDBQT format."""
    obabel_bin = shutil.which("obabel")
    if obabel_bin:
        cmd = [obabel_bin, "-ipdb", str(clean_pdb), "-opdbqt", "-O", str(out_pdbqt), "-xr", "--partialcharge", "gasteiger"]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode == 0 and out_pdbqt.exists() and out_pdbqt.stat().st_size > 1000:
            return "OpenBabel Rigid Receptor (-xr, Gasteiger)"

    out_lines: List[str] = []
    for ln in clean_pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if ln.startswith(("ATOM  ", "HETATM")):
            atom_name = ln[12:16].strip()
            elem = (ln[76:78].strip() if len(ln) >= 78 else "").upper()
            if not elem:
                elem = "".join([c for c in atom_name if c.isalpha()])[:1].upper() or "C"
            if elem == "H":
                if not any(tag in atom_name.upper() for tag in ("HN", "HO", "HS", "H1", "H2", "H3", "HG", "HH", "HD", "HE", "HZ")):
                    continue
                ad4_type = "HD"
            else:
                ad4_type = AD4_PROTON_MAP.get(elem, "C")
            base = ln[:54].ljust(54)
            out_lines.append(f"{base}  1.00  0.00    +0.000 {ad4_type:<2}")
        elif ln.startswith("TER"):
            out_lines.append("TER")
    out_pdbqt.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    return "AnCore AD4 Receptor Typer"


def prepare_reference_ligand_pdbqt(crystal_lig_pdb: Path, work_dir: Path) -> Tuple[Path, Path, str]:
    """
    Prepare co-crystal ligand SDF and FLEXIBLE PDBQT (with active rotatable bond torsions)
    while preserving 3D crystal coordinates.
    """
    sdf_path = work_dir / "crystal_ligand.sdf"
    pdbqt_path = work_dir / "crystal_ligand.pdbqt"

    obabel_bin = shutil.which("obabel")
    if obabel_bin:
        subprocess.run(
            [obabel_bin, "-ipdb", str(crystal_lig_pdb), "-osdf", "-O", str(sdf_path), "-h"],
            capture_output=True, text=True,
        )
        res = subprocess.run(
            [obabel_bin, "-ipdb", str(crystal_lig_pdb), "-opdbqt", "-O", str(pdbqt_path), "-h", "--partialcharge", "gasteiger"],
            capture_output=True, text=True,
        )
        if res.returncode == 0 and pdbqt_path.exists() and pdbqt_path.stat().st_size > 100:
            return sdf_path, pdbqt_path, "OpenBabel Flexible Torsion Tree (Gasteiger)"

    try:
        from rdkit import Chem
        from meeko import MoleculePreparation, PDBQTWriterLegacy

        mol = Chem.MolFromPDBFile(str(crystal_lig_pdb), removeHs=False, sanitize=True)
        if mol is not None:
            mol_h = Chem.AddHs(mol, addCoords=True)
            w = Chem.SDWriter(str(sdf_path))
            w.write(mol_h)
            w.close()

            prep = MoleculePreparation()
            setups = prep.prepare(mol_h)
            for setup in setups:
                pdbqt_str, is_ok, _ = PDBQTWriterLegacy.write_string(setup)
                if is_ok:
                    pdbqt_path.write_text(pdbqt_str, encoding="utf-8")
                    return sdf_path, pdbqt_path, "RDKit + Meeko Flexible Torsion Tree"
    except Exception:
        pass

    out_lines = ["ROOT"]
    for ln in crystal_lig_pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if ln.startswith(("ATOM  ", "HETATM")):
            elem = (ln[76:78].strip() if len(ln) >= 78 else ln[12:16].strip()[:1]).upper() or "C"
            ad4_type = AD4_PROTON_MAP.get(elem, "C")
            base = ln[:54].ljust(54)
            out_lines.append(f"{base}  1.00  0.00    +0.000 {ad4_type:<2}")
    out_lines.extend(["ENDROOT", "TORSDOF 0"])
    pdbqt_path.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    return sdf_path, pdbqt_path, "Native Fallback"


def detect_pocket_residues(
    clean_receptor_pdb: Path,
    lig_coords: List[Tuple[float, float, float]],
    cutoff: float = 5.0,
) -> Tuple[List[Dict[str, Any]], List[int]]:
    """Identify all receptor residues with any heavy atom within `cutoff` A of the crystal ligand."""
    res_min_dist: Dict[Tuple[str, int, str], float] = {}
    cutoff_sq = cutoff * cutoff

    for ln in clean_receptor_pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if not ln.startswith("ATOM  "):
            continue
        elem = (ln[76:78].strip() if len(ln) >= 78 else ln[12:16].strip()[:1]).upper()
        if elem == "H":
            continue
        try:
            rx = float(ln[30:38])
            ry = float(ln[38:46])
            rz = float(ln[46:54])
            resname = ln[17:20].strip()
            chain = ln[21:22].strip() or "A"
            resnum = int(ln[22:26].strip())
        except ValueError:
            continue

        key = (chain, resnum, resname)
        for lx, ly, lz in lig_coords:
            d2 = (rx - lx) ** 2 + (ry - ly) ** 2 + (rz - lz) ** 2
            if d2 <= cutoff_sq:
                d = math.sqrt(d2)
                if key not in res_min_dist or d < res_min_dist[key]:
                    res_min_dist[key] = d

    hydrophobic = {"ALA", "VAL", "LEU", "ILE", "MET", "PRO", "GLY"}
    aromatic = {"PHE", "TYR", "TRP", "HIS", "HID", "HIE", "HIP"}
    polar = {"SER", "THR", "ASN", "GLN", "CYS", "CYX"}

    pocket_records: List[Dict[str, Any]] = []
    for (chain, resnum, resname), dist in sorted(res_min_dist.items(), key=lambda x: x[1]):
        if resname in aromatic:
            chem_class = "Aromatic (pi-Stacking / Hydrophobic)"
        elif resname in hydrophobic:
            chem_class = "Aliphatic Hydrophobic"
        elif resname in polar:
            chem_class = "Polar Uncharged (H-Bond)"
        else:
            chem_class = "Charged / Salt-Bridge"

        pocket_records.append({
            "chain": chain,
            "res_num": resnum,
            "res_name": resname,
            "label": f"{resname}{resnum}",
            "min_dist_A": round(dist, 2),
            "chem_class": chem_class,
        })

    pocket_res_ids = sorted({r["res_num"] for r in pocket_records})
    return pocket_records, pocket_res_ids


def compute_search_box(
    lig_coords: List[Tuple[float, float, float]],
    padding: float = 4.0,
    min_size: float = 16.0,
    max_size: float = 24.0,
) -> Dict[str, float]:
    """Compute 3D bounding box centered on the co-crystallized ligand."""
    xs = [c[0] for c in lig_coords]
    ys = [c[1] for c in lig_coords]
    zs = [c[2] for c in lig_coords]

    cx = (min(xs) + max(xs)) / 2.0
    cy = (min(ys) + max(ys)) / 2.0
    cz = (min(zs) + max(zs)) / 2.0

    sx = max(min_size, min(max_size, (max(xs) - min(xs)) + 2.0 * padding))
    sy = max(min_size, min(max_size, (max(ys) - min(ys)) + 2.0 * padding))
    sz = max(min_size, min(max_size, (max(zs) - min(zs)) + 2.0 * padding))

    return {
        "center_x": round(cx, 3),
        "center_y": round(cy, 3),
        "center_z": round(cz, 3),
        "size_x": round(sx, 2),
        "size_y": round(sy, 2),
        "size_z": round(sz, 2),
        "padding_A": round(padding, 2),
        "volume_A3": round(sx * sy * sz, 1),
    }


# ==============================================================================
# SINGLE-TARGET & MULTI-BENCHMARK SUITE RUNNERS
# ==============================================================================
def _curate_single_target(
    proj_path: Path,
    pdb_code: str,
    chain_id: Optional[str] = None,
    ref_ligand: Optional[str] = None,
    ph: float = 7.4,
    box_padding: float = 4.0,
    pocket_cutoff: float = 5.0,
    progress_cb: Optional[Callable[[int, str], None]] = None,
) -> Dict[str, Any]:
    """Curate a single PDB target and save its manifest."""
    t0 = time.time()
    pdb_code = pdb_code.strip().upper()
    bench_info = BENCHMARK_TARGETS.get(pdb_code, {})

    chosen_chain = (chain_id or bench_info.get("chain") or "A").strip().upper()
    chosen_lig = (ref_ligand or bench_info.get("ligand") or "").strip().upper() or None

    target_dir = proj_path / "benchmarks" / pdb_code
    target_dir.mkdir(parents=True, exist_ok=True)

    if progress_cb:
        progress_cb(15, f"[{pdb_code}] Downloading crystal structure from RCSB PDB...")
    raw_pdb = download_pdb(pdb_code, target_dir)

    if progress_cb:
        progress_cb(35, f"[{pdb_code}] Isolating Chain {chosen_chain} and crystal ligand...")
    isolated_pdb, crystal_lig_pdb, det_chain, det_lig, lig_coords = isolate_receptor_and_ligand(
        raw_pdb, target_dir, target_chain=chosen_chain, ref_ligand=chosen_lig
    )

    if progress_cb:
        progress_cb(55, f"[{pdb_code}] Repairing heavy atoms & protonating at pH {ph:.1f}...")
    clean_pdb = target_dir / "receptor_clean.pdb"
    repair_engine = repair_and_protonate_receptor(isolated_pdb, clean_pdb, ph=ph)

    if progress_cb:
        progress_cb(75, f"[{pdb_code}] Parameterizing receptor & flexible ligand PDBQTs...")
    receptor_pdbqt = target_dir / "receptor.pdbqt"
    rec_pdbqt_engine = receptor_pdb_to_pdbqt(clean_pdb, receptor_pdbqt)
    crystal_lig_sdf, crystal_lig_pdbqt, lig_pdbqt_engine = prepare_reference_ligand_pdbqt(crystal_lig_pdb, target_dir)

    if progress_cb:
        progress_cb(90, f"[{pdb_code}] Mapping {pocket_cutoff:.1f} A pocket & 3D search box...")
    pocket_residues, pocket_res_ids = detect_pocket_residues(clean_pdb, lig_coords, cutoff=pocket_cutoff)
    box = compute_search_box(lig_coords, padding=box_padding)

    vina_box_txt = target_dir / "vina_box.txt"
    vina_box_txt.write_text(
        f"receptor = {receptor_pdbqt}\n"
        f"center_x = {box['center_x']}\n"
        f"center_y = {box['center_y']}\n"
        f"center_z = {box['center_z']}\n"
        f"size_x = {box['size_x']}\n"
        f"size_y = {box['size_y']}\n"
        f"size_z = {box['size_z']}\n",
        encoding="utf-8",
    )

    for src_file in (clean_pdb, receptor_pdbqt, crystal_lig_pdb, crystal_lig_pdbqt, vina_box_txt):
        if src_file.exists():
            shutil.copyfile(src_file, proj_path / src_file.name)

    elapsed = round(time.time() - t0, 2)
    manifest: Dict[str, Any] = {
        "stage": 1,
        "stage_name": "Target Curation, Repair & 3D Search Box",
        "status": "COMPLETED",
        "elapsed_seconds": elapsed,
        "project_dir": str(proj_path),
        "target_dir": str(target_dir),
        "pdb_id": pdb_code,
        "target_name": bench_info.get("name", f"PDB {pdb_code}"),
        "ligand_name": bench_info.get("ligand_name", det_lig),
        "resolution": bench_info.get("resolution", "Crystal"),
        "benchmark_source": bench_info.get("benchmark", "Custom RCSB Target"),
        "pocket_architecture": bench_info.get("pocket_type", "Crystallographic Binding Pocket"),
        "chain_id": det_chain,
        "ligand_id": det_lig,
        "ph": ph,
        "engines": {
            "receptor_repair": repair_engine,
            "receptor_pdbqt": rec_pdbqt_engine,
            "ligand_pdbqt": lig_pdbqt_engine,
        },
        "box": box,
        "pocket_cutoff_A": pocket_cutoff,
        "num_pocket_residues": len(pocket_residues),
        "pocket_residue_ids": pocket_res_ids,
        "pocket_residues": pocket_residues,
        "files": {
            "raw_pdb": str(raw_pdb),
            "receptor_clean_pdb": str(clean_pdb),
            "receptor_pdb": str(clean_pdb),
            "receptor_pdbqt": str(receptor_pdbqt),
            "crystal_ligand_pdb": str(crystal_lig_pdb),
            "crystal_ligand_sdf": str(crystal_lig_sdf),
            "crystal_ligand_pdbqt": str(crystal_lig_pdbqt),
            "vina_box_txt": str(vina_box_txt),
        },
    }

    (target_dir / "stage1_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (proj_path / "stage1_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def run_stage1(
    project_dir: str = "./",
    pdb_id: str = "ALL",
    chain_id: Optional[str] = None,
    ref_ligand: Optional[str] = None,
    ligand_id: Optional[str] = None,
    ph: float = 7.4,
    box_padding: float = 4.0,
    pocket_cutoff: float = 5.0,
    run_all_benchmarks: bool = False,
    progress_cb: Optional[Callable[[int, str], None]] = None,
    **kwargs: Any,
) -> Dict[str, Any]:
    """
    Execute Stage 1 Target Curation.
    If `pdb_id="ALL"` (default) or `run_all_benchmarks=True`, curates all 5 gold-standard
    benchmark proteins (1IEP, 3ERT, 1UYG, 1M17, 1STP) and prints a 100% completion summary.
    """
    t_start = time.time()
    proj_path = Path(project_dir).expanduser().resolve()
    proj_path.mkdir(parents=True, exist_ok=True)

    req_code = pdb_id.strip().upper()
    if req_code == "ALL" or run_all_benchmarks:
        codes = list(BENCHMARK_TARGETS.keys())
        suite_results: Dict[str, Dict[str, Any]] = {}
        n_total = len(codes)
        for idx, code in enumerate(codes):
            base_pct = int((idx / n_total) * 100)
            span = int(100 / n_total)

            def _sub_cb(p: int, m: str, b=base_pct, s=span) -> None:
                overall = min(99, b + int((p / 100.0) * s))
                logger.info(f"[{overall}%] {m}")
                if progress_cb:
                    progress_cb(overall, m)

            m_data = _curate_single_target(
                proj_path=proj_path,
                pdb_code=code,
                ph=ph,
                box_padding=box_padding,
                pocket_cutoff=pocket_cutoff,
                progress_cb=_sub_cb,
            )
            suite_results[code] = m_data

        suite_summary = {
            "stage": 1,
            "suite_completed": True,
            "targets": suite_results,
            "active_pdb_id": "1IEP",
        }
        (proj_path / "stage1_suite_manifest.json").write_text(json.dumps(suite_summary, indent=2), encoding="utf-8")

        total_time = round(time.time() - t_start, 2)
        done_msg = f"[100%] Stage 1 Complete: All {n_total} Benchmark Targets (1IEP, 3ERT, 1UYG, 1M17, 1STP) Curated in {total_time}s."
        logger.info(done_msg)
        print("=" * 80)
        print(f" [100%] ANCORE STAGE 1 BENCHMARK SUITE CURATED SUCCESSFULLY ({total_time}s)")
        print("=" * 80)
        for code in codes:
            m = suite_results.get(code, {})
            b = m.get("box", {})
            print(
                f" * {code:<4} | {m.get('target_name', ''):<26} | Ligand: {m.get('ligand_id', ''):<4} | "
                f"Pocket: {m.get('num_pocket_residues', 0):>2} res | "
                f"Box: {b.get('size_x', 0):>4.1f}x{b.get('size_y', 0):>4.1f}x{b.get('size_z', 0):>4.1f} A"
            )
        print("=" * 80)
        if progress_cb:
            progress_cb(100, done_msg)
        return suite_results["1IEP"]
    else:
        chosen_lig = ref_ligand or ligand_id
        manifest = _curate_single_target(
            proj_path=proj_path,
            pdb_code=req_code,
            chain_id=chain_id,
            ref_ligand=chosen_lig,
            ph=ph,
            box_padding=box_padding,
            pocket_cutoff=pocket_cutoff,
            progress_cb=progress_cb,
        )
        total_time = round(time.time() - t_start, 2)
        box = manifest["box"]
        done_msg = (
            f"[100%] Stage 1 Complete | Target: {req_code} ({manifest['num_pocket_residues']} pocket residues) | "
            f"Box: ({box['size_x']} x {box['size_y']} x {box['size_z']} A) in {total_time}s."
        )
        logger.info(done_msg)
        if progress_cb:
            progress_cb(100, done_msg)
        return manifest


def load_all_stage1_manifests(proj_path: Path) -> Dict[str, Dict[str, Any]]:
    """Scan `proj_path/benchmarks/<PDB>/stage1_manifest.json` and return all prepared targets."""
    found: Dict[str, Dict[str, Any]] = {}
    for code in BENCHMARK_TARGETS:
        m_file = proj_path / "benchmarks" / code / "stage1_manifest.json"
        if m_file.exists():
            try:
                found[code] = json.loads(m_file.read_text(encoding="utf-8"))
            except Exception:
                pass
    root_m = proj_path / "stage1_manifest.json"
    if root_m.exists():
        try:
            d = json.loads(root_m.read_text(encoding="utf-8"))
            c = d.get("pdb_id", "").upper()
            if c and c not in found:
                found[c] = d
        except Exception:
            pass
    return found


# ==============================================================================
# ANCORE STUDIO 3Dmol.js VIEWPORT (Stage 1)
# ==============================================================================
def build_stage1_3d_html(manifest: Dict[str, Any], highlight_res: Optional[int] = None) -> str:
    """Build AnCore Dark Studio 3Dmol.js viewport with HUD overlay, cavity surface toggle, and box."""
    files = manifest.get("files", {})
    rec_path = Path(files.get("receptor_clean_pdb") or files.get("receptor_pdb", ""))
    lig_path = Path(files.get("crystal_ligand_pdb", ""))

    rec_pdb_str = rec_path.read_text(encoding="utf-8", errors="replace") if rec_path.exists() else ""
    lig_pdb_str = lig_path.read_text(encoding="utf-8", errors="replace") if lig_path.exists() else ""
    pocket_ids = manifest.get("pocket_residue_ids", [])
    box = manifest.get("box", {})
    pdb_id = manifest.get("pdb_id", "1IEP")
    target_name = manifest.get("target_name", pdb_id)
    lig_id = manifest.get("ligand_id", "LIG")
    lig_name = manifest.get("ligand_name", lig_id)
    bench_src = manifest.get("benchmark_source", "Benchmark")
    pocket_arch = manifest.get("pocket_architecture", "")

    return f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
  * {{ box-sizing: border-box; margin: 0; padding: 0; font-family: "Helvetica Neue", Helvetica, Arial, sans-serif; }}
  html, body {{ width: 100%; height: 100%; background: #090D14; overflow: hidden; user-select: none; }}
  #viewport {{ width: 100%; height: 100%; position: relative; }}
  .hud-card {{
    position: absolute; top: 12px; left: 14px; z-index: 20;
    background: rgba(15, 19, 26, 0.86); border: 1px solid #262E3D;
    border-left: 3px solid #38BDF8; padding: 8px 12px; border-radius: 6px;
    color: #F1F5F9; font-size: 11px; line-height: 1.45; max-width: 420px;
    pointer-events: none;
  }}
  .hud-title {{ font-size: 13px; font-weight: 800; color: #38BDF8; letter-spacing: 0.3px; }}
  .hud-sub {{ color: #94A3B8; font-size: 10.5px; }}
  .toolbar {{
    position: absolute; top: 12px; right: 14px; z-index: 20;
    display: flex; gap: 5px; background: rgba(15, 19, 26, 0.90);
    padding: 5px 7px; border-radius: 6px; border: 1px solid #262E3D;
  }}
  .toolbar button {{
    background: #1A1F2B; color: #CBD5E1; border: 1px solid #334155;
    padding: 5px 10px; border-radius: 4px; font-size: 11px; font-weight: 600;
    cursor: pointer; transition: all 0.12s;
  }}
  .toolbar button.active {{ background: #0284C7; color: #FFFFFF; border-color: #38BDF8; }}
  .toolbar button:hover {{ background: #262E3D; color: #FFFFFF; }}
  .legend {{
    position: absolute; bottom: 12px; left: 14px; z-index: 20;
    background: rgba(15, 19, 26, 0.88); border: 1px solid #262E3D;
    padding: 7px 12px; border-radius: 6px; color: #CBD5E1; font-size: 11px;
    display: flex; gap: 14px; align-items: center; pointer-events: none;
  }}
  .dot {{ display: inline-block; width: 9px; height: 9px; border-radius: 50%; margin-right: 5px; vertical-align: middle; }}
</style>
</head>
<body>
<div id="viewport"></div>
<div class="hud-card">
  <div class="hud-title">{pdb_id} — {target_name}</div>
  <div><b>Co-Crystal Ligand:</b> {lig_id} ({lig_name}) &nbsp;|&nbsp; <b>Set:</b> {bench_src}</div>
  <div class="hud-sub">{pocket_arch}</div>
</div>
<div class="toolbar">
  <button id="btnRec" class="active" onclick="toggleRec()">Receptor Ribbon</button>
  <button id="btnPock" class="active" onclick="togglePock()">5A Pocket ({len(pocket_ids)})</button>
  <button id="btnSurf" onclick="toggleSurf()">Pocket Surface</button>
  <button id="btnLig" class="active" onclick="toggleLig()">Crystal ({lig_id})</button>
  <button id="btnBox" class="active" onclick="toggleBox()">3D Grid Box</button>
  <button onclick="resetCam()">Reset View</button>
</div>
<div class="legend">
  <span><i class="dot" style="background:#64748B;"></i>Receptor ({pdb_id})</span>
  <span><i class="dot" style="background:#F59E0B;"></i>5.0 A Pocket Residues</span>
  <span><i class="dot" style="background:#10B981;"></i>Crystal Ligand ({lig_id})</span>
  <span><i class="dot" style="background:#38BDF8;"></i>Grid ({box.get('size_x', 0)} x {box.get('size_y', 0)} x {box.get('size_z', 0)} A)</span>
</div>
<script>
  const recData = {json.dumps(rec_pdb_str)};
  const ligData = {json.dumps(lig_pdb_str)};
  const pocketIds = {json.dumps(pocket_ids)};
  const box = {json.dumps(box)};
  const highlightRes = {json.dumps(highlight_res)};

  let viewer = null;
  let recModel = null;
  let ligModel = null;
  let surfObj = null;
  let boxShapes = [];
  let vis = {{ rec: true, pock: true, surf: false, lig: true, box: true }};

  function applyStyles() {{
    if (!viewer) return;
    if (recModel) {{
      recModel.setStyle({{}}, {{}});
      if (vis.rec) {{
        recModel.setStyle({{}}, {{ cartoon: {{ color: "#475569", opacity: 0.72 }} }});
      }}
      if (vis.pock && pocketIds.length > 0) {{
        recModel.setStyle({{ resi: pocketIds }}, {{
          cartoon: {{ color: vis.rec ? "#F59E0B" : "transparent", opacity: 0.88 }},
          stick: {{ colorscheme: "orangeCarbon", radius: 0.14 }}
        }});
      }}
      if (highlightRes !== null) {{
        recModel.setStyle({{ resi: [highlightRes] }}, {{
          stick: {{ colorscheme: "magentaCarbon", radius: 0.28 }}
        }});
      }}
    }}
    if (ligModel) {{
      ligModel.setStyle({{}}, vis.lig ? {{ stick: {{ colorscheme: "greenCarbon", radius: 0.28 }} }} : {{}});
    }}
    if (surfObj !== null) {{
      viewer.removeSurface(surfObj);
      surfObj = null;
    }}
    if (vis.surf && recModel && pocketIds.length > 0) {{
      surfObj = viewer.addSurface($3Dmol.SurfaceType.VDW, {{
        opacity: 0.38, color: "#0EA5E9"
      }}, {{ model: recModel, resi: pocketIds }});
    }}
    boxShapes.forEach(s => viewer.removeShape(s));
    boxShapes = [];
    if (vis.box && box.center_x !== undefined) {{
      const s1 = viewer.addBox({{
        center: {{ x: box.center_x, y: box.center_y, z: box.center_z }},
        dimensions: {{ w: box.size_x, h: box.size_y, d: box.size_z }},
        color: "#38BDF8", opacity: 0.12, wireframe: false
      }});
      const s2 = viewer.addBox({{
        center: {{ x: box.center_x, y: box.center_y, z: box.center_z }},
        dimensions: {{ w: box.size_x, h: box.size_y, d: box.size_z }},
        color: "#38BDF8", opacity: 0.90, wireframe: true
      }});
      boxShapes.push(s1, s2);
    }}
    viewer.render();
  }}

  function setBtn(id, state) {{
    const el = document.getElementById(id);
    if (el) el.className = state ? "active" : "";
  }}
  function toggleRec() {{ vis.rec = !vis.rec; setBtn("btnRec", vis.rec); applyStyles(); }}
  function togglePock() {{ vis.pock = !vis.pock; setBtn("btnPock", vis.pock); applyStyles(); }}
  function toggleSurf() {{ vis.surf = !vis.surf; setBtn("btnSurf", vis.surf); applyStyles(); }}
  function toggleLig() {{ vis.lig = !vis.lig; setBtn("btnLig", vis.lig); applyStyles(); }}
  function toggleBox() {{ vis.box = !vis.box; setBtn("btnBox", vis.box); applyStyles(); }}
  function resetCam() {{
    if (!viewer) return;
    if (highlightRes !== null && recModel) {{
      viewer.zoomTo({{ model: recModel, resi: [highlightRes] }}, 400);
    }} else if (ligModel) {{
      viewer.zoomTo({{ model: ligModel }}, 400);
      viewer.zoom(0.76);
    }} else {{
      viewer.zoomTo();
    }}
    viewer.render();
  }}

  window.addEventListener("load", function() {{
    viewer = $3Dmol.createViewer("viewport", {{ backgroundColor: "#090D14" }});
    if (recData) recModel = viewer.addModel(recData, "pdb");
    if (ligData) ligModel = viewer.addModel(ligData, "pdb");
    applyStyles();
    resetCam();
  }});
</script>
</body>
</html>"""


# ==============================================================================
# ANCORE STUDIO NATIVE macOS GUI (Stage 1)
# ==============================================================================
_STAGE1_WINDOW_REF = None


def launch_stage1_gui(default_project_dir: str = "./", default_pdb: str = "1IEP") -> Any:
    """Launch the AnCore Dark Studio GUI for Stage 1."""
    global _STAGE1_WINDOW_REF
    try:
        from PySide6.QtCore import QObject, Qt, QThread, QUrl, Signal
        from PySide6.QtGui import QColor, QDoubleValidator
        from PySide6.QtWidgets import (
            QApplication, QFrame, QGridLayout, QHBoxLayout, QHeaderView,
            QLabel, QLineEdit, QMainWindow, QMessageBox, QProgressBar,
            QPushButton, QSplitter, QTabWidget, QTableWidget,
            QTableWidgetItem, QVBoxLayout, QWidget,
        )
        try:
            from PySide6.QtWebEngineWidgets import QWebEngineView
            HAS_WEBENGINE = True
        except ImportError:
            HAS_WEBENGINE = False
    except ImportError as exc:
        raise RuntimeError("PySide6 is required for the macOS GUI. Install with: pip install PySide6") from exc

    class Stage1Worker(QObject):
        progress = Signal(int, str)
        finished = Signal(dict)
        failed = Signal(str)

        def __init__(self, params: Dict[str, Any]) -> None:
            super().__init__()
            self.params = params

        def run(self) -> None:
            try:
                res = run_stage1(
                    **self.params,
                    progress_cb=lambda pct, msg: self.progress.emit(pct, msg),
                )
                self.finished.emit(res)
            except Exception as exc:
                self.failed.emit(str(exc))

    class AnCoreStage1Window(QMainWindow):
        def __init__(self) -> None:
            super().__init__()
            self.setWindowTitle("AnCore Studio — Stage 1: Multi-Benchmark Target Curation & 3D Grid Generation")
            self.resize(1440, 900)
            self.proj_path = Path(default_project_dir).expanduser().resolve()
            self.active_pdb = default_pdb.strip().upper() if default_pdb.strip().upper() in BENCHMARK_TARGETS else "1IEP"
            self.manifests_by_pdb: Dict[str, Dict[str, Any]] = {}
            self._worker_thread: Optional[QThread] = None

            self.setStyleSheet("""
                QMainWindow, QWidget {
                    background-color: #141820;
                    color: #F1F5F9;
                    font-family: "Helvetica Neue", Helvetica, Arial, sans-serif;
                    font-size: 12px;
                }
                QFrame#cardFrame {
                    background-color: #1A1F2B;
                    border: 1px solid #283040;
                    border-radius: 7px;
                }
                QLabel#cardTitle {
                    color: #38BDF8;
                    font-size: 11px;
                    font-weight: 800;
                    letter-spacing: 0.6px;
                }
                QLabel#fieldLabel {
                    color: #94A3B8;
                    font-size: 10.5px;
                    font-weight: 700;
                }
                QLineEdit {
                    background-color: #10141C;
                    color: #F8FAFC;
                    border: 1px solid #334155;
                    border-radius: 5px;
                    padding: 6px 8px;
                    font-size: 12px;
                    font-weight: 600;
                }
                QLineEdit:focus {
                    border: 1px solid #38BDF8;
                    background-color: #0D1118;
                }
                QTableWidget {
                    background-color: #141820;
                    alternate-background-color: #181D27;
                    color: #E2E8F0;
                    gridline-color: #262E3D;
                    border: 1px solid #283040;
                    border-radius: 5px;
                    font-size: 11.5px;
                    selection-background-color: #0369A1;
                    selection-color: #FFFFFF;
                }
                QHeaderView::section {
                    background-color: #10141C;
                    color: #94A3B8;
                    font-weight: 700;
                    font-size: 11px;
                    border: none;
                    border-bottom: 1px solid #283040;
                    border-right: 1px solid #1E2430;
                    padding: 6px 8px;
                }
                QTabWidget::pane {
                    border: 1px solid #283040;
                    background-color: #1A1F2B;
                    border-radius: 5px;
                }
                QTabBar::tab {
                    background-color: #10141C;
                    color: #94A3B8;
                    border: 1px solid #283040;
                    border-bottom: none;
                    padding: 7px 16px;
                    margin-right: 4px;
                    border-top-left-radius: 5px;
                    border-top-right-radius: 5px;
                    font-weight: 700;
                    font-size: 11.5px;
                }
                QTabBar::tab:selected {
                    background-color: #1A1F2B;
                    color: #38BDF8;
                    border-top: 2px solid #38BDF8;
                }
                QProgressBar {
                    background-color: #10141C;
                    border: 1px solid #283040;
                    border-radius: 4px;
                    text-align: center;
                    color: #F8FAFC;
                    font-size: 10px;
                    font-weight: 700;
                }
                QProgressBar::chunk {
                    background-color: #0284C7;
                    border-radius: 3px;
                }
                QSplitter::handle {
                    background-color: #262E3D;
                }
            """)

            central = QWidget()
            self.setCentralWidget(central)
            root_layout = QVBoxLayout(central)
            root_layout.setContentsMargins(0, 0, 0, 0)
            root_layout.setSpacing(0)

            # 1. TOP COMMAND RIBBON (52px)
            top_bar = QFrame()
            top_bar.setFixedHeight(52)
            top_bar.setStyleSheet("background-color: #0D1117; border-bottom: 1px solid #262E3D;")
            top_layout = QHBoxLayout(top_bar)
            top_layout.setContentsMargins(18, 0, 18, 0)

            brand = QLabel(
                "<span style='color:#38BDF8;font-size:18px;font-weight:900;letter-spacing:0.5px;'>AnCore</span> "
                "<span style='color:#64748B;font-size:12px;font-weight:700;'>SCIENTIFIC STUDIO</span>"
                "<span style='color:#334155;font-size:14px;'> &nbsp;|&nbsp; </span>"
                "<span style='color:#E2E8F0;font-size:13px;font-weight:700;'>"
                "STAGE 1: STRUCTURE CURATION, PROTONATION &amp; 3D SEARCH GRID</span>"
            )
            top_layout.addWidget(brand)
            top_layout.addStretch()

            self.suite_badge = QLabel("BENCHMARK SUITE: 0 / 5 PREPARED")
            self.suite_badge.setStyleSheet(
                "background-color: #1A1F2B; color: #38BDF8; font-weight: 800; font-size: 11.5px; "
                "padding: 6px 14px; border-radius: 5px; border: 1px solid #0284C7;"
            )
            top_layout.addWidget(self.suite_badge)
            root_layout.addWidget(top_bar)

            # 2. MAIN VERTICAL SPLITTER
            v_splitter = QSplitter(Qt.Orientation.Vertical)
            v_splitter.setHandleWidth(3)
            root_layout.addWidget(v_splitter, 1)

            h_splitter = QSplitter(Qt.Orientation.Horizontal)
            h_splitter.setHandleWidth(3)
            v_splitter.addWidget(h_splitter)

            # --- LEFT STUDIO CONTROL DOCK ---
            left_dock = QWidget()
            left_dock.setMinimumWidth(375)
            left_dock.setMaximumWidth(450)
            left_layout = QVBoxLayout(left_dock)
            left_layout.setContentsMargins(10, 10, 10, 10)
            left_layout.setSpacing(10)

            # CARD 1: BENCHMARK SUITE NAVIGATOR (All 5 Targets Visible!)
            nav_card = QFrame()
            nav_card.setObjectName("cardFrame")
            nav_layout = QVBoxLayout(nav_card)
            nav_layout.setContentsMargins(10, 10, 10, 10)
            nav_layout.setSpacing(6)

            nav_title = QLabel("BENCHMARK SUITE NAVIGATOR (CLICK ANY TARGET)")
            nav_title.setObjectName("cardTitle")
            nav_layout.addWidget(nav_title)

            self.nav_table = QTableWidget(len(BENCHMARK_TARGETS), 4)
            self.nav_table.setHorizontalHeaderLabels(["PDB", "Target Receptor", "Ligand", "Status"])
            self.nav_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
            self.nav_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
            self.nav_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
            self.nav_table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)
            self.nav_table.verticalHeader().setVisible(False)
            self.nav_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
            self.nav_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
            self.nav_table.setAlternatingRowColors(True)
            self.nav_table.setFixedHeight(172)
            self.nav_table.cellClicked.connect(self.on_nav_target_clicked)
            nav_layout.addWidget(self.nav_table)
            left_layout.addWidget(nav_card)

            # CARD 2: CLEAN 2-COLUMN PARAMETER GRID
            param_card = QFrame()
            param_card.setObjectName("cardFrame")
            p_layout = QVBoxLayout(param_card)
            p_layout.setContentsMargins(10, 10, 10, 10)
            p_layout.setSpacing(8)

            p_title = QLabel("TARGET CURATION & GRID PARAMETERS")
            p_title.setObjectName("cardTitle")
            p_layout.addWidget(p_title)

            lbl_proj = QLabel("WORKSPACE PROJECT DIRECTORY")
            lbl_proj.setObjectName("fieldLabel")
            self.proj_input = QLineEdit(str(self.proj_path))
            p_layout.addWidget(lbl_proj)
            p_layout.addWidget(self.proj_input)

            grid = QGridLayout()
            grid.setHorizontalSpacing(10)
            grid.setVerticalSpacing(6)

            l_pdb = QLabel("ACTIVE PDB ID"); l_pdb.setObjectName("fieldLabel")
            l_chain = QLabel("TARGET CHAIN"); l_chain.setObjectName("fieldLabel")
            l_lig = QLabel("CO-CRYSTAL LIGAND"); l_lig.setObjectName("fieldLabel")
            l_ph = QLabel("PROTONATION pH"); l_ph.setObjectName("fieldLabel")
            l_pad = QLabel("BOX PADDING (A)"); l_pad.setObjectName("fieldLabel")
            l_cut = QLabel("POCKET CUTOFF (A)"); l_cut.setObjectName("fieldLabel")

            self.pdb_input = QLineEdit(self.active_pdb)
            self.chain_input = QLineEdit("A")
            self.lig_input = QLineEdit("STI")
            self.ph_input = QLineEdit("7.4")
            self.ph_input.setValidator(QDoubleValidator(4.0, 10.0, 1))
            self.pad_input = QLineEdit("4.0")
            self.pad_input.setValidator(QDoubleValidator(1.5, 12.0, 1))
            self.cut_input = QLineEdit("5.0")
            self.cut_input.setValidator(QDoubleValidator(3.0, 10.0, 1))

            grid.addWidget(l_pdb, 0, 0); grid.addWidget(l_chain, 0, 1)
            grid.addWidget(self.pdb_input, 1, 0); grid.addWidget(self.chain_input, 1, 1)
            grid.addWidget(l_lig, 2, 0); grid.addWidget(l_ph, 2, 1)
            grid.addWidget(self.lig_input, 3, 0); grid.addWidget(self.ph_input, 3, 1)
            grid.addWidget(l_pad, 4, 0); grid.addWidget(l_cut, 4, 1)
            grid.addWidget(self.pad_input, 5, 0); grid.addWidget(self.cut_input, 5, 1)

            p_layout.addLayout(grid)
            left_layout.addWidget(param_card)

            # CARD 3: EXECUTION CONTROLS & KPI SUMMARY
            exec_card = QFrame()
            exec_card.setObjectName("cardFrame")
            e_layout = QVBoxLayout(exec_card)
            e_layout.setContentsMargins(10, 10, 10, 10)
            e_layout.setSpacing(8)

            self.btn_run_all = QPushButton("Prepare All 5 Benchmark Targets (Full Suite)")
            self.btn_run_all.setStyleSheet(
                "QPushButton { background-color: #0284C7; color: #FFFFFF; font-weight: 800; "
                "padding: 10px; border-radius: 5px; font-size: 12.5px; border: 1px solid #38BDF8; }"
                "QPushButton:hover { background-color: #0369A1; }"
                "QPushButton:disabled { background-color: #334155; color: #64748B; border-color: #334155; }"
            )
            self.btn_run_all.clicked.connect(lambda: self.start_curation(run_all=True))
            e_layout.addWidget(self.btn_run_all)

            self.btn_run_single = QPushButton("Prepare Selected Target Only")
            self.btn_run_single.setStyleSheet(
                "QPushButton { background-color: #262E3D; color: #E2E8F0; font-weight: 700; "
                "padding: 8px; border-radius: 5px; font-size: 12px; border: 1px solid #475569; }"
                "QPushButton:hover { background-color: #334155; color: #FFFFFF; }"
                "QPushButton:disabled { background-color: #1E2430; color: #64748B; }"
            )
            self.btn_run_single.clicked.connect(lambda: self.start_curation(run_all=False))
            e_layout.addWidget(self.btn_run_single)

            self.progress_bar = QProgressBar()
            self.progress_bar.setRange(0, 100)
            self.progress_bar.setValue(0)
            self.progress_bar.setFixedHeight(14)
            e_layout.addWidget(self.progress_bar)

            self.status_lbl = QLabel("Ready. Click 'Prepare All 5 Benchmark Targets' to curate the benchmark suite.")
            self.status_lbl.setWordWrap(True)
            self.status_lbl.setStyleSheet("color: #94A3B8; font-size: 11px;")
            e_layout.addWidget(self.status_lbl)

            kpi_row = QHBoxLayout()
            kpi_row.setSpacing(8)
            self.kpi_pocket = QLabel("<b>5A POCKET</b><br><span style='font-size:15px;color:#38BDF8;'>—</span>")
            self.kpi_pocket.setAlignment(Qt.AlignmentFlag.AlignCenter)
            self.kpi_pocket.setStyleSheet("background:#10141C;border:1px solid #283040;border-radius:5px;padding:7px;")
            self.kpi_volume = QLabel("<b>GRID VOLUME</b><br><span style='font-size:15px;color:#10B981;'>—</span>")
            self.kpi_volume.setAlignment(Qt.AlignmentFlag.AlignCenter)
            self.kpi_volume.setStyleSheet("background:#10141C;border:1px solid #283040;border-radius:5px;padding:7px;")
            kpi_row.addWidget(self.kpi_pocket)
            kpi_row.addWidget(self.kpi_volume)
            e_layout.addLayout(kpi_row)

            left_layout.addWidget(exec_card)
            left_layout.addStretch()
            h_splitter.addWidget(left_dock)

            # --- RIGHT WIDESCREEN 3D MOLECULAR VIEWPORT ---
            if HAS_WEBENGINE:
                self.web_view = QWebEngineView()
                self.web_view.setHtml(
                    "<body style='background:#090D14;color:#64748B;display:flex;align-items:center;"
                    "justify-content:center;height:100vh;font-family:Helvetica,Arial,sans-serif;'>"
                    "<div>Click <b>Prepare All 5 Benchmark Targets</b> to curate and visualize the benchmark suite.</div></body>"
                )
                h_splitter.addWidget(self.web_view)
            else:
                fb = QLabel("Install PySide6-Addons for interactive 3Dmol.js WebGL viewport.")
                fb.setAlignment(Qt.AlignmentFlag.AlignCenter)
                h_splitter.addWidget(fb)

            h_splitter.setSizes([390, 1050])

            # --- BOTTOM PROJECT TABLE (TABBED) ---
            self.bottom_tabs = QTabWidget()
            self.bottom_tabs.setContentsMargins(8, 4, 8, 8)

            self.suite_table = QTableWidget(len(BENCHMARK_TARGETS), 9)
            self.suite_table.setHorizontalHeaderLabels([
                "PDB ID", "Target Receptor", "Benchmark Set", "Chain", "Ligand",
                "5A Pocket Residues", "Grid Center (x, y, z)", "Grid Size (A)", "Curation Status"
            ])
            self.suite_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
            self.suite_table.verticalHeader().setVisible(False)
            self.suite_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
            self.suite_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
            self.suite_table.setAlternatingRowColors(True)
            self.suite_table.cellClicked.connect(self.on_suite_table_clicked)
            self.bottom_tabs.addTab(self.suite_table, "Benchmark Suite Overview (All 5 Targets)")

            self.pocket_table = QTableWidget(0, 5)
            self.pocket_table.setHorizontalHeaderLabels([
                "Residue", "Chain", "Residue #", "Min Contact Distance (A)", "Chemical Interaction Class"
            ])
            self.pocket_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
            self.pocket_table.verticalHeader().setVisible(False)
            self.pocket_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
            self.pocket_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
            self.pocket_table.setAlternatingRowColors(True)
            self.pocket_table.cellClicked.connect(self.on_pocket_residue_clicked)
            self.bottom_tabs.addTab(self.pocket_table, "Active Target 5.0 A Binding Pocket Residues")

            v_splitter.addWidget(self.bottom_tabs)
            v_splitter.setSizes([610, 238])

            self.refresh_all_views()

        def refresh_all_views(self) -> None:
            self.proj_path = Path(self.proj_input.text().strip() or "./").expanduser().resolve()
            self.manifests_by_pdb = load_all_stage1_manifests(self.proj_path)

            codes = list(BENCHMARK_TARGETS.keys())
            n_ready = sum(1 for c in codes if c in self.manifests_by_pdb)
            self.suite_badge.setText(f"BENCHMARK SUITE: {n_ready} / {len(codes)} TARGETS PREPARED")
            if n_ready == len(codes):
                self.suite_badge.setStyleSheet(
                    "background-color: #064E3B; color: #34D399; font-weight: 800; font-size: 11.5px; "
                    "padding: 6px 14px; border-radius: 5px; border: 1px solid #10B981;"
                )

            for row_idx, code in enumerate(codes):
                meta = BENCHMARK_TARGETS[code]
                m = self.manifests_by_pdb.get(code)
                is_ready = m is not None

                nav_items = [
                    QTableWidgetItem(code),
                    QTableWidgetItem(meta["name"]),
                    QTableWidgetItem(meta["ligand"]),
                    QTableWidgetItem("READY" if is_ready else "PENDING"),
                ]
                for c_i, it in enumerate(nav_items):
                    if c_i in (0, 2, 3):
                        it.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    if c_i == 0:
                        it.setForeground(QColor("#38BDF8"))
                    if c_i == 3:
                        it.setForeground(QColor("#10B981" if is_ready else "#F59E0B"))
                    self.nav_table.setItem(row_idx, c_i, it)

                if is_ready and m:
                    box = m.get("box", {})
                    center_str = f"({box.get('center_x')}, {box.get('center_y')}, {box.get('center_z')})"
                    size_str = f"{box.get('size_x')} x {box.get('size_y')} x {box.get('size_z')} A"
                    n_res_str = f"{m.get('num_pocket_residues', 0)} residues"
                    stat_str = f"CURATED ({m.get('elapsed_seconds', 0)}s)"
                else:
                    center_str = "—"
                    size_str = "—"
                    n_res_str = "—"
                    stat_str = "Not Curated Yet (Click Prepare All 5)"

                suite_items = [
                    QTableWidgetItem(code),
                    QTableWidgetItem(meta["name"]),
                    QTableWidgetItem(meta["benchmark"]),
                    QTableWidgetItem(meta["chain"]),
                    QTableWidgetItem(f"{meta['ligand']} ({meta['ligand_name']})"),
                    QTableWidgetItem(n_res_str),
                    QTableWidgetItem(center_str),
                    QTableWidgetItem(size_str),
                    QTableWidgetItem(stat_str),
                ]
                for c_i, it in enumerate(suite_items):
                    it.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    if c_i == 0:
                        it.setForeground(QColor("#38BDF8"))
                    if c_i == 8:
                        it.setForeground(QColor("#10B981" if is_ready else "#94A3B8"))
                    self.suite_table.setItem(row_idx, c_i, it)

            self.select_active_target(self.active_pdb)

        def select_active_target(self, pdb_code: str) -> None:
            pdb_code = pdb_code.strip().upper()
            self.active_pdb = pdb_code
            meta = BENCHMARK_TARGETS.get(pdb_code, {})
            self.pdb_input.setText(pdb_code)
            if meta:
                self.chain_input.setText(meta["chain"])
                self.lig_input.setText(meta["ligand"])

            codes = list(BENCHMARK_TARGETS.keys())
            if pdb_code in codes:
                r_idx = codes.index(pdb_code)
                self.nav_table.selectRow(r_idx)
                self.suite_table.selectRow(r_idx)

            m = self.manifests_by_pdb.get(pdb_code)
            if m:
                self.populate_active_manifest(m)

        def on_nav_target_clicked(self, row: int, _col: int) -> None:
            codes = list(BENCHMARK_TARGETS.keys())
            if 0 <= row < len(codes):
                self.select_active_target(codes[row])

        def on_suite_table_clicked(self, row: int, _col: int) -> None:
            codes = list(BENCHMARK_TARGETS.keys())
            if 0 <= row < len(codes):
                self.select_active_target(codes[row])

        def on_pocket_residue_clicked(self, row: int, _col: int) -> None:
            m = self.manifests_by_pdb.get(self.active_pdb)
            if not m or not HAS_WEBENGINE:
                return
            residues = m.get("pocket_residues", [])
            if 0 <= row < len(residues):
                res_num = residues[row].get("res_num")
                html_str = build_stage1_3d_html(m, highlight_res=res_num)
                self.web_view.setHtml(html_str, QUrl("https://3Dmol.org/"))

        def start_curation(self, run_all: bool = True) -> None:
            self.btn_run_all.setEnabled(False)
            self.btn_run_single.setEnabled(False)
            self.progress_bar.setValue(5)

            target_pdb = "ALL" if run_all else (self.pdb_input.text().strip().upper() or "1IEP")
            if not run_all:
                self.active_pdb = target_pdb

            params = {
                "project_dir": self.proj_input.text().strip() or "./",
                "pdb_id": target_pdb,
                "chain_id": self.chain_input.text().strip().upper() or "A",
                "ref_ligand": self.lig_input.text().strip().upper() or None,
                "ph": float(self.ph_input.text() or 7.4),
                "box_padding": float(self.pad_input.text() or 4.0),
                "pocket_cutoff": float(self.cut_input.text() or 5.0),
                "run_all_benchmarks": run_all,
            }

            self._worker_thread = QThread()
            self._worker = Stage1Worker(params)
            self._worker.moveToThread(self._worker_thread)
            self._worker_thread.started.connect(self._worker.run)
            self._worker.progress.connect(self.on_progress)
            self._worker.finished.connect(self.on_finished)
            self._worker.failed.connect(self.on_failed)
            self._worker_thread.start()

        def on_progress(self, pct: int, msg: str) -> None:
            self.progress_bar.setValue(pct)
            self.status_lbl.setText(msg)

        def on_finished(self, _manifest: Dict[str, Any]) -> None:
            if self._worker_thread:
                self._worker_thread.quit()
                self._worker_thread.wait()
            self.btn_run_all.setEnabled(True)
            self.btn_run_single.setEnabled(True)
            self.refresh_all_views()

        def on_failed(self, err_msg: str) -> None:
            if self._worker_thread:
                self._worker_thread.quit()
                self._worker_thread.wait()
            self.btn_run_all.setEnabled(True)
            self.btn_run_single.setEnabled(True)
            self.status_lbl.setText(f"Error: {err_msg}")
            QMessageBox.critical(self, "Stage 1 Error", err_msg)

        def populate_active_manifest(self, manifest: Dict[str, Any]) -> None:
            pdb_id = manifest.get("pdb_id", "1IEP")
            n_res = manifest.get("num_pocket_residues", 0)
            box = manifest.get("box", {})
            vol = box.get("volume_A3", 0.0)

            self.kpi_pocket.setText(
                f"<span style='font-size:10px;color:#94A3B8;font-weight:700;'>{pdb_id} 5A POCKET</span><br>"
                f"<span style='font-size:15px;font-weight:800;color:#38BDF8;'>{n_res} Residues</span>"
            )
            self.kpi_volume.setText(
                f"<span style='font-size:10px;color:#94A3B8;font-weight:700;'>GRID VOLUME</span><br>"
                f"<span style='font-size:15px;font-weight:800;color:#10B981;'>{vol:,.0f} A³</span>"
            )

            self.bottom_tabs.setTabText(1, f"Active Target ({pdb_id}) — 5.0 A Binding Pocket Residues ({n_res})")

            residues = manifest.get("pocket_residues", [])
            self.pocket_table.setRowCount(len(residues))
            for r_idx, r in enumerate(residues):
                dist_val = r.get("min_dist_A", 9.9)
                items = [
                    QTableWidgetItem(r.get("label", "")),
                    QTableWidgetItem(str(r.get("chain", "A"))),
                    QTableWidgetItem(str(r.get("res_num", ""))),
                    QTableWidgetItem(f"{dist_val:.2f} A"),
                    QTableWidgetItem(r.get("chem_class", "")),
                ]
                for c_idx, it in enumerate(items):
                    it.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    if dist_val <= 3.3 and c_idx in (0, 3):
                        it.setForeground(QColor("#38BDF8"))
                    self.pocket_table.setItem(r_idx, c_idx, it)

            if HAS_WEBENGINE:
                html_str = build_stage1_3d_html(manifest)
                self.web_view.setHtml(html_str, QUrl("https://3Dmol.org/"))

    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    win = AnCoreStage1Window()
    _STAGE1_WINDOW_REF = win
    win.show()
    win.raise_()
    win.activateWindow()
    logger.info("Stage 1 AnCore Studio GUI opened. Close the window when finished to return to Jupyter.")
    app.exec()
    return win


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="AnCore Stage 1 — Multi-Benchmark Target Curation")
    parser.add_argument("--project-dir", type=str, default="./", help="Project directory")
    parser.add_argument("--pdb-id", type=str, default="ALL", help="RCSB PDB ID or 'ALL' for all 5 benchmarks")
    parser.add_argument("--ph", type=float, default=7.4, help="Protonation pH")
    parser.add_argument("--padding", type=float, default=4.0, help="Box padding in Angstroms")
    parser.add_argument("--gui", action="store_true", help="Launch native macOS GUI")

    clean_argv = [a for a in sys.argv[1:] if not a.startswith("-f") and "kernel-" not in a and not a.endswith(".json")]
    args, _ = parser.parse_known_args(clean_argv)

    if args.gui:
        launch_stage1_gui(args.project_dir, args.pdb_id if args.pdb_id != "ALL" else "1IEP")
    else:
        run_stage1(project_dir=args.project_dir, pdb_id=args.pdb_id, ph=args.ph, box_padding=args.padding)
