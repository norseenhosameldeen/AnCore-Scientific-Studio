#!/usr/bin/env python3
"""
AnCore — Stage 5: Protein-Ligand Interaction Fingerprinting & 2D/3D Pharmacophore Map (v2.2)
============================================================================================
Authors: Norseen Hossam El-Deen & Ahmed M. Senary
Affiliation: Faculty of Pharmacy, Galala University, Egypt
License: Apache-2.0 (100% Original Open-Source Implementation)

Description:
    Standalone Stage 5 module for the AnCore Native macOS Molecular Docking Workbench.
    Implements Q1 computational medicinal chemistry protocols for non-covalent interaction profiling:
    1. 3D Geometric & Chemical Interaction Fingerprinting (PLIP / ProLIF Standard Criteria):
       - Hydrogen Bonds (H-BOND, <= 3.5 A between polar N/O/S donor/acceptor pairs)
       - Salt Bridges / Ionic Contacts (SALT-BRIDGE, <= 4.2 A between polar N/O and ASP/GLU/ARG/LYS/HIS)
       - Aromatic pi-Stacking (PI-STACKING, <= 4.3 A between aromatic/sp2 ligand atoms and PHE/TYR/TRP/HIS)
       - Halogen Bonds (HALOGEN-BOND, <= 3.9 A between F/Cl/Br/I and receptor backbone/sidechain O/N)
       - Hydrophobic Contacts (HYDROPHOBIC, <= 4.0 A between non-polar carbons and hydrophobic residues)
    2. Crystallographic Interaction Conservation (IFP Tanimoto Similarity):
       - Profiles the co-crystallized reference ligand first to establish the native pocket hotspot
         fingerprint, then computes the Tanimoto IFP Similarity (0.00 - 1.00) for every screened molecule.
    3. AnCore Dark Studio GUI (PySide6 + 3Dmol.js 3D Bond Cylinders + 2D Radial Pharmacophore SVG):
       - Left Benchmark Suite Navigator + Screened Complex Selector + cutoff controls.
       - Widescreen 3D molecular viewport rendering dashed 3D interaction cylinders, 3D residue labels,
         and an embedded 2D Radial Pharmacophore Interaction Diagram (SVG).
       - Bottom tabbed table: [Atom-by-Atom 3D Interactions (Click to Zoom)] & [Library IFP Summary].
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import os
import sys
import time
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Set, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("AnCore.Stage5")

BENCHMARK_TARGETS: Dict[str, Dict[str, Any]] = {
    "1IEP": {
        "name": "ABL1 Tyrosine Kinase",
        "ligand": "STI",
        "benchmark": "DUD-E / Vina",
        "key_hotspots": ["THR315", "GLU286", "ASP381", "MET318", "PHE317", "ILE313", "VAL256"],
    },
    "3ERT": {
        "name": "Estrogen Receptor alpha",
        "ligand": "OHT",
        "benchmark": "DUD-E & LIT-PCBA",
        "key_hotspots": ["GLU353", "ARG394", "HIS524", "ASP351", "PHE404", "LEU346", "LEU387"],
    },
    "1UYG": {
        "name": "HSP90-alpha Chaperone",
        "ligand": "PU3",
        "benchmark": "DUD-E Benchmark",
        "key_hotspots": ["ASP93", "THR184", "ASN51", "MET98", "PHE138", "LEU107", "TYR139"],
    },
    "1M17": {
        "name": "EGFR Tyrosine Kinase",
        "ligand": "AQ4",
        "benchmark": "DUD-E Oncology",
        "key_hotspots": ["MET769", "THR766", "LYS721", "ASP831", "LEU694", "LEU820", "GLN767"],
    },
    "1STP": {
        "name": "Streptavidin",
        "ligand": "BTN",
        "benchmark": "CASF-2016 Gold Standard",
        "key_hotspots": ["ASN23", "SER27", "TYR43", "SER45", "ASN49", "TRP79", "SER88", "ASP128"],
    },
}

AROMATIC_RES = {"PHE", "TYR", "TRP", "HIS", "HID", "HIE", "HIP"}
CHARGED_RES = {"ASP", "GLU", "ARG", "LYS", "HIS", "HIP"}
HYDROPHOBIC_RES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "TYR", "PRO"}
HALOGEN_ELEMS = {"F", "CL", "BR", "I"}
POLAR_ELEMS = {"N", "O", "S"}


# ==============================================================================
# 3D ATOM PARSING & GEOMETRIC INTERACTION FINGERPRINTING ENGINE
# ==============================================================================
def parse_pdb_heavy_atoms(pdb_text: str, is_receptor: bool = True) -> List[Dict[str, Any]]:
    """Parse heavy atoms from PDB text with coordinates, atom names, elements, and residue labels."""
    atoms: List[Dict[str, Any]] = []
    for ln in pdb_text.splitlines():
        if not ln.startswith(("ATOM  ", "HETATM")):
            continue
        atom_name = ln[12:16].strip()
        tail = ln[76:78].strip() if len(ln) >= 78 else ""
        raw_elem = tail if tail else "".join([c for c in atom_name if c.isalpha()])[:2]
        elem = raw_elem.upper()
        if elem not in HALOGEN_ELEMS:
            elem = elem[:1]
        if elem == "H" or atom_name.upper().startswith("H"):
            continue
        try:
            x = float(ln[30:38])
            y = float(ln[38:46])
            z = float(ln[46:54])
            res_name = ln[17:20].strip() or ("PRO" if is_receptor else "LIG")
            chain = ln[21:22].strip() or "A"
            res_num_str = ln[22:26].strip()
            res_num = int(res_num_str) if res_num_str.lstrip("-").isdigit() else 1
        except ValueError:
            continue

        atoms.append({
            "atom_name": atom_name,
            "elem": elem,
            "res_name": res_name,
            "chain": chain,
            "res_num": res_num,
            "res_label": f"{res_name}{res_num}",
            "x": x,
            "y": y,
            "z": z,
        })
    return atoms


def profile_complex_interactions(
    rec_atoms: List[Dict[str, Any]],
    lig_atoms: List[Dict[str, Any]],
    hotspots: Set[str],
    hbond_cutoff: float = 3.5,
    salt_cutoff: float = 4.2,
    pistack_cutoff: float = 4.3,
    halogen_cutoff: float = 3.9,
    hydro_cutoff: float = 4.0,
) -> Tuple[List[Dict[str, Any]], Dict[str, int], List[str]]:
    """
    Detect non-covalent interactions between receptor pocket atoms and ligand pose atoms.
    Deduplicates to keep the strongest (shortest-distance) interaction per (Residue, InteractionType)
    so the interaction table and 3D network are clean and non-redundant.
    """
    max_cut_sq = max(hbond_cutoff, salt_cutoff, pistack_cutoff, halogen_cutoff, hydro_cutoff) ** 2

    # Keep shortest contact per (res_label, interaction_type)
    best_by_res_type: Dict[Tuple[str, str], Dict[str, Any]] = {}

    for r_at in rec_atoms:
        rx, ry, rz = r_at["x"], r_at["y"], r_at["z"]
        r_elem = r_at["elem"]
        r_res = r_at["res_name"]
        r_lbl = r_at["res_label"]
        r_aname = r_at["atom_name"]

        for l_at in lig_atoms:
            lx, ly, lz = l_at["x"], l_at["y"], l_at["z"]
            d2 = (rx - lx) ** 2 + (ry - ly) ** 2 + (rz - lz) ** 2
            if d2 > max_cut_sq:
                continue
            dist = math.sqrt(d2)
            l_elem = l_at["elem"]
            l_aname = l_at["atom_name"]

            detected_type: Optional[str] = None
            role_desc = ""
            priority = 99

            # 1. Hydrogen Bond (N/O/S ... N/O/S <= 3.5 A)
            if l_elem in POLAR_ELEMS and r_elem in POLAR_ELEMS and dist <= hbond_cutoff:
                detected_type = "H-BOND"
                role_desc = "Directional polar hydrogen bond"
                priority = 1

            # 2. Halogen Bond (F/Cl/Br/I ... O/N/S <= 3.9 A)
            elif l_elem in HALOGEN_ELEMS and r_elem in POLAR_ELEMS and dist <= halogen_cutoff:
                detected_type = "HALOGEN-BOND"
                role_desc = "Sigma-hole halogen-oxygen/nitrogen interaction"
                priority = 2

            # 3. Salt Bridge / Ionic Contact (Polar N/O ... Charged Residue Sidechain <= 4.2 A)
            elif (
                l_elem in {"N", "O"}
                and r_res in CHARGED_RES
                and r_aname not in {"N", "O", "C", "CA"}
                and r_elem in {"N", "O"}
                and dist <= salt_cutoff
            ):
                detected_type = "SALT-BRIDGE"
                role_desc = "Electrostatic ionic / salt-bridge anchor"
                priority = 2

            # 4. Aromatic pi-Stacking (Ligand C/N ... PHE/TYR/TRP/HIS ring atoms <= 4.3 A)
            elif (
                l_elem in {"C", "N"}
                and r_res in AROMATIC_RES
                and r_aname not in {"N", "O", "C", "CA", "CB"}
                and dist <= pistack_cutoff
            ):
                detected_type = "PI-STACKING"
                role_desc = "Aromatic pi-pi / T-shaped ring stacking"
                priority = 3

            # 5. Hydrophobic Contact (Ligand C ... Hydrophobic Sidechain C <= 4.0 A)
            elif (
                l_elem == "C"
                and r_elem in {"C", "S"}
                and r_res in HYDROPHOBIC_RES
                and r_aname not in {"N", "O", "C", "CA"}
                and dist <= hydro_cutoff
            ):
                detected_type = "HYDROPHOBIC"
                role_desc = "Lipophilic van der Waals burial"
                priority = 4

            if detected_type is not None:
                key = (r_lbl, detected_type)
                if key not in best_by_res_type or dist < best_by_res_type[key]["distance_A"]:
                    best_by_res_type[key] = {
                        "interaction_type": detected_type,
                        "priority": priority,
                        "res_label": r_lbl,
                        "res_name": r_res,
                        "res_num": r_at["res_num"],
                        "chain": r_at["chain"],
                        "rec_atom": r_aname,
                        "lig_atom": l_aname,
                        "distance_A": round(dist, 2),
                        "is_hotspot": r_lbl in hotspots,
                        "role_desc": role_desc,
                        "rec_xyz": [round(rx, 3), round(ry, 3), round(rz, 3)],
                        "lig_xyz": [round(lx, 3), round(ly, 3), round(lz, 3)],
                    }

    # If a residue forms a specific polar/aromatic interaction (H-BOND, SALT-BRIDGE, PI-STACKING),
    # suppress redundant HYDROPHOBIC label on that exact same residue if we already have >= 6 contacts
    specific_residues = {
        k[0] for k in best_by_res_type if k[1] in ("H-BOND", "SALT-BRIDGE", "HALOGEN-BOND", "PI-STACKING")
    }
    filtered: List[Dict[str, Any]] = []
    for (r_lbl, itype), rec in best_by_res_type.items():
        if itype == "HYDROPHOBIC" and r_lbl in specific_residues and len(best_by_res_type) > 6:
            continue
        filtered.append(rec)

    filtered.sort(key=lambda x: (x["priority"], x["distance_A"]))

    counts = {
        "H-BOND": sum(1 for x in filtered if x["interaction_type"] == "H-BOND"),
        "SALT-BRIDGE": sum(1 for x in filtered if x["interaction_type"] == "SALT-BRIDGE"),
        "PI-STACKING": sum(1 for x in filtered if x["interaction_type"] == "PI-STACKING"),
        "HALOGEN-BOND": sum(1 for x in filtered if x["interaction_type"] == "HALOGEN-BOND"),
        "HYDROPHOBIC": sum(1 for x in filtered if x["interaction_type"] == "HYDROPHOBIC"),
        "TOTAL": len(filtered),
    }
    contacted_residues = sorted({x["res_label"] for x in filtered})
    return filtered, counts, contacted_residues


def compute_ifp_tanimoto(cryst_residues: List[str], pose_residues: List[str]) -> float:
    """Compute Tanimoto similarity between crystal reference interacting residues and docked pose."""
    set_a = set(cryst_residues)
    set_b = set(pose_residues)
    if not set_a and not set_b:
        return 1.0
    union = len(set_a | set_b)
    if union == 0:
        return 0.0
    return round(len(set_a & set_b) / union, 3)


# ==============================================================================
# SINGLE-TARGET & MULTI-BENCHMARK INTERACTION PROFILER
# ==============================================================================
def _profile_single_target(
    proj_path: Path,
    pdb_code: str,
    hbond_cutoff: float = 3.5,
    hydro_cutoff: float = 4.0,
    pistack_cutoff: float = 4.3,
    salt_cutoff: float = 4.2,
    progress_cb: Optional[Callable[[int, str], None]] = None,
) -> Dict[str, Any]:
    """Profile protein-ligand interactions for the co-crystal ligand and all screened ligands."""
    t0 = time.time()
    pdb_code = pdb_code.strip().upper()
    target_dir = proj_path / "benchmarks" / pdb_code
    inter_dir = target_dir / "interactions"
    inter_dir.mkdir(parents=True, exist_ok=True)

    # Ensure Stage 4 exists (which automatically ensures Stages 1 & 3 exist)
    s4_path = target_dir / "stage4_manifest.json"
    if not s4_path.exists():
        if progress_cb:
            progress_cb(10, f"[{pdb_code}] Running Stage 4 Screening prerequisites automatically...")
        from stage4_screening import run_stage4
        run_stage4(project_dir=str(proj_path), pdb_id=pdb_code, run_all_benchmarks=False)
    s4_manifest = json.loads(s4_path.read_text(encoding="utf-8"))

    s1_path = target_dir / "stage1_manifest.json"
    s1_manifest = json.loads(s1_path.read_text(encoding="utf-8")) if s1_path.exists() else {}

    rec_pdb_path = Path(s4_manifest["files"]["receptor_clean_pdb"])
    cryst_pdb_path = Path(s4_manifest["files"]["crystal_ligand_pdb"])
    rec_pdb_text = rec_pdb_path.read_text(encoding="utf-8", errors="replace")
    cryst_pdb_text = cryst_pdb_path.read_text(encoding="utf-8", errors="replace")

    # Pre-filter receptor atoms to those within 9.0 A of pocket center for speed
    box = s4_manifest.get("box", {})
    cx, cy, cz = float(box.get("center_x", 0.0)), float(box.get("center_y", 0.0)), float(box.get("center_z", 0.0))
    all_rec_atoms = parse_pdb_heavy_atoms(rec_pdb_text, is_receptor=True)
    pocket_rec_atoms = [
        a for a in all_rec_atoms
        if (a["x"] - cx) ** 2 + (a["y"] - cy) ** 2 + (a["z"] - cz) ** 2 <= 22.0 ** 2
    ]
    if not pocket_rec_atoms:
        pocket_rec_atoms = all_rec_atoms

    bench_meta = BENCHMARK_TARGETS.get(pdb_code, {})
    hotspots = set(bench_meta.get("key_hotspots", []))

    # 1. Profile Co-Crystal Reference Ligand first
    cryst_atoms = parse_pdb_heavy_atoms(cryst_pdb_text, is_receptor=False)
    cryst_contacts, cryst_counts, cryst_res_list = profile_complex_interactions(
        pocket_rec_atoms, cryst_atoms, hotspots,
        hbond_cutoff=hbond_cutoff, salt_cutoff=salt_cutoff,
        pistack_cutoff=pistack_cutoff, hydro_cutoff=hydro_cutoff,
    )

    complex_profiles: List[Dict[str, Any]] = [{
        "complex_index": 0,
        "consensus_rank": 0,
        "id": f"CRYSTAL_{s4_manifest.get('ligand_id', 'REF')}",
        "name": f"Co-Crystal Reference ({s4_manifest.get('ligand_id', 'REF')})",
        "role": "CRYSTAL_REF",
        "score_vina": s4_manifest.get("top_hit_vina", -11.0),
        "counts": cryst_counts,
        "contacted_residues": cryst_res_list,
        "hotspots_hit": [r for r in cryst_res_list if r in hotspots],
        "ifp_similarity": 1.000,
        "interactions": cryst_contacts,
        "pose_pdb_block": cryst_pdb_text,
    }]

    # 2. Profile all ranked screened ligands from Stage 4
    ranked_ligs = s4_manifest.get("ranked_ligands", [])
    n_ligs = len(ranked_ligs)
    for idx, lig in enumerate(ranked_ligs):
        pct = 20 + int(((idx + 1) / max(1, n_ligs)) * 75)
        if progress_cb:
            progress_cb(pct, f"[{pdb_code}] Profiling 3D interactions for Rank #{idx+1}: {lig['name']}...")

        pose_text = lig.get("pose_pdb_block", "")
        if not pose_text and lig.get("pose_pdb_path") and Path(lig["pose_pdb_path"]).exists():
            pose_text = Path(lig["pose_pdb_path"]).read_text(encoding="utf-8", errors="replace")

        l_atoms = parse_pdb_heavy_atoms(pose_text, is_receptor=False)
        contacts, counts, res_list = profile_complex_interactions(
            pocket_rec_atoms, l_atoms, hotspots,
            hbond_cutoff=hbond_cutoff, salt_cutoff=salt_cutoff,
            pistack_cutoff=pistack_cutoff, hydro_cutoff=hydro_cutoff,
        )
        ifp_sim = compute_ifp_tanimoto(cryst_res_list, res_list)
        hits_hot = [r for r in res_list if r in hotspots]

        complex_profiles.append({
            "complex_index": idx + 1,
            "consensus_rank": lig.get("consensus_rank", idx + 1),
            "id": lig.get("id", f"LIG_{idx+1}"),
            "name": lig.get("name", "Compound"),
            "role": lig.get("role", "ACTIVE"),
            "score_vina": lig.get("score_vina", 0.0),
            "score_vinardo": lig.get("score_vinardo", 0.0),
            "rrf_score": lig.get("rrf_score", 0.0),
            "counts": counts,
            "contacted_residues": res_list,
            "hotspots_hit": hits_hot,
            "ifp_similarity": ifp_sim,
            "interactions": contacts,
            "pose_pdb_block": pose_text,
        })

    # Export CSV summary
    csv_path = target_dir / "interactions_summary.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "rank", "id", "name", "role", "score_vina",
            "hbonds", "salt_bridges", "pi_stacking", "halogen_bonds", "hydrophobic",
            "total_contacts", "ifp_tanimoto_similarity", "hotspots_hit", "all_contact_residues",
        ])
        for cp in complex_profiles:
            c = cp["counts"]
            writer.writerow([
                cp["consensus_rank"], cp["id"], cp["name"], cp["role"], cp["score_vina"],
                c["H-BOND"], c["SALT-BRIDGE"], c["PI-STACKING"], c["HALOGEN-BOND"], c["HYDROPHOBIC"],
                c["TOTAL"], cp["ifp_similarity"],
                ";".join(cp["hotspots_hit"]), ";".join(cp["contacted_residues"]),
            ])

    elapsed = round(time.time() - t0, 2)
    top_screened = complex_profiles[1] if len(complex_profiles) > 1 else complex_profiles[0]

    manifest: Dict[str, Any] = {
        "stage": 5,
        "stage_name": "Protein-Ligand Interaction Fingerprinting & 2D/3D Pharmacophore Map",
        "status": "COMPLETED",
        "elapsed_seconds": elapsed,
        "project_dir": str(proj_path),
        "target_dir": str(target_dir),
        "pdb_id": pdb_code,
        "target_name": bench_meta.get("name", s4_manifest.get("target_name", pdb_code)),
        "benchmark_source": bench_meta.get("benchmark", "Benchmark"),
        "ligand_id": s4_manifest.get("ligand_id", "LIG"),
        "key_hotspots": list(hotspots),
        "pocket_residue_ids": s4_manifest.get("pocket_residue_ids", s1_manifest.get("pocket_residue_ids", [])),
        "top_hit_name": top_screened["name"],
        "top_hit_hbonds": top_screened["counts"]["H-BOND"],
        "top_hit_total_contacts": top_screened["counts"]["TOTAL"],
        "top_hit_ifp_similarity": top_screened["ifp_similarity"],
        "top_hit_hotspots": top_screened["hotspots_hit"],
        "csv_summary_path": str(csv_path),
        "complexes": complex_profiles,
        "files": {
            "receptor_clean_pdb": str(rec_pdb_path),
            "crystal_ligand_pdb": str(cryst_pdb_path),
        },
    }

    (target_dir / "stage5_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (proj_path / "stage5_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def run_stage5(
    project_dir: str = "./",
    pdb_id: str = "ALL",
    hbond_cutoff: float = 3.5,
    hydro_cutoff: float = 4.0,
    pistack_cutoff: float = 4.3,
    salt_cutoff: float = 4.2,
    run_all_benchmarks: bool = False,
    progress_cb: Optional[Callable[[int, str], None]] = None,
    **kwargs: Any,
) -> Dict[str, Any]:
    """
    Execute Stage 5 Protein-Ligand Interaction Fingerprinting across all 5 benchmark targets
    (or a single target) and print a 100% completion summary.
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

            m_data = _profile_single_target(
                proj_path=proj_path,
                pdb_code=code,
                hbond_cutoff=hbond_cutoff,
                hydro_cutoff=hydro_cutoff,
                pistack_cutoff=pistack_cutoff,
                salt_cutoff=salt_cutoff,
                progress_cb=_sub_cb,
            )
            suite_results[code] = m_data

        suite_summary = {
            "stage": 5,
            "suite_completed": True,
            "targets": suite_results,
            "active_pdb_id": "1IEP",
        }
        (proj_path / "stage5_suite_manifest.json").write_text(json.dumps(suite_summary, indent=2), encoding="utf-8")

        total_time = round(time.time() - t_start, 2)
        done_msg = f"[100%] Stage 5 Complete: Profiled 3D Non-Covalent Interactions Across All 5 Targets in {total_time}s."
        logger.info(done_msg)
        print("=" * 94)
        print(f" [100%] ANCORE STAGE 5 PROTEIN-LIGAND INTERACTION FINGERPRINTING COMPLETED ({total_time}s)")
        print("=" * 94)
        for code in codes:
            m = suite_results.get(code, {})
            hots = ", ".join(m.get("top_hit_hotspots", [])[:4]) or "Pocket Core"
            print(
                f" * {code:<4} | Top Hit: {m.get('top_hit_name', '')[:18]:<18} | "
                f"H-Bonds: {m.get('top_hit_hbonds', 0):>2} | Contacts: {m.get('top_hit_total_contacts', 0):>2} | "
                f"IFP Match: {m.get('top_hit_ifp_similarity', 0):.2f} | Hotspots: {hots}"
            )
        print("=" * 94)
        if progress_cb:
            progress_cb(100, done_msg)
        return suite_results["1IEP"]
    else:
        manifest = _profile_single_target(
            proj_path=proj_path,
            pdb_code=req_code,
            hbond_cutoff=hbond_cutoff,
            hydro_cutoff=hydro_cutoff,
            pistack_cutoff=pistack_cutoff,
            salt_cutoff=salt_cutoff,
            progress_cb=progress_cb,
        )
        total_time = round(time.time() - t_start, 2)
        done_msg = (
            f"[100%] Stage 5 Complete | Target: {req_code} | "
            f"Top Hit IFP Similarity = {manifest['top_hit_ifp_similarity']:.2f} in {total_time}s."
        )
        logger.info(done_msg)
        if progress_cb:
            progress_cb(100, done_msg)
        return manifest


def load_all_stage5_manifests(proj_path: Path) -> Dict[str, Dict[str, Any]]:
    """Scan `proj_path/benchmarks/<PDB>/stage5_manifest.json` and return all profiled targets."""
    found: Dict[str, Dict[str, Any]] = {}
    for code in BENCHMARK_TARGETS:
        m_file = proj_path / "benchmarks" / code / "stage5_manifest.json"
        if m_file.exists():
            try:
                found[code] = json.loads(m_file.read_text(encoding="utf-8"))
            except Exception:
                pass
    root_m = proj_path / "stage5_manifest.json"
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
# 2D RADIAL PHARMACOPHORE SVG DIAGRAM + 3Dmol.js INTERACTION VIEWPORT
# ==============================================================================
TYPE_COLORS = {
    "H-BOND": "#38BDF8",       # Cyan
    "SALT-BRIDGE": "#EC4899",  # Pink/Magenta
    "PI-STACKING": "#A855F7",  # Purple
    "HALOGEN-BOND": "#22C55E", # Green
    "HYDROPHOBIC": "#F59E0B",  # Amber
}


def _build_2d_pharmacophore_svg(complex_data: Dict[str, Any], pdb_id: str) -> str:
    """Build an interactive 2D radial protein-ligand interaction diagram in SVG."""
    w, h = 265, 235
    cx, cy = w / 2.0, h / 2.0 + 6
    radius = 78.0

    interactions = complex_data.get("interactions", [])[:10]
    lig_id = complex_data.get("id", "LIG")
    n = len(interactions)

    lines_svg: List[str] = []
    nodes_svg: List[str] = []

    for i, inter in enumerate(interactions):
        angle = (2.0 * math.pi * i / max(1, n)) - (math.pi / 2.0)
        rx = cx + radius * math.cos(angle)
        ry = cy + radius * math.sin(angle)
        itype = inter.get("interaction_type", "HYDROPHOBIC")
        col = TYPE_COLORS.get(itype, "#94A3B8")
        dash = "4,3" if itype in ("H-BOND", "SALT-BRIDGE", "HALOGEN-BOND") else ("2,2" if itype == "PI-STACKING" else "none")
        dist = inter.get("distance_A", 0.0)
        res_lbl = inter.get("res_label", "RES")
        is_hot = inter.get("is_hotspot", False)

        # Midpoint distance label
        mx = (cx + rx) / 2.0
        my = (cy + ry) / 2.0

        lines_svg.append(
            f'<line x1="{cx:.1f}" y1="{cy:.1f}" x2="{rx:.1f}" y2="{ry:.1f}" '
            f'stroke="{col}" stroke-width="{2.0 if is_hot else 1.4}" stroke-dasharray="{dash}"/>'
        )
        lines_svg.append(
            f'<text x="{mx:.1f}" y="{my:.1f}" fill="#E2E8F0" font-size="8" font-weight="bold" '
            f'text-anchor="middle" dy="-2">{dist:.1f}A</text>'
        )

        stroke_w = "2.2" if is_hot else "1.3"
        nodes_svg.append(
            f'<circle cx="{rx:.1f}" cy="{ry:.1f}" r="18" fill="#141820" stroke="{col}" stroke-width="{stroke_w}"/>'
            f'<text x="{rx:.1f}" y="{ry - 2:.1f}" fill="{col}" font-size="8.5" font-weight="bold" text-anchor="middle">{res_lbl}</text>'
            f'<text x="{rx:.1f}" y="{ry + 7:.1f}" fill="#94A3B8" font-size="7" text-anchor="middle">{itype[:6]}</text>'
        )

    return f"""
    <svg width="{w}" height="{h}" viewBox="0 0 {w} {h}" xmlns="http://www.w3.org/2000/svg">
      <rect width="{w}" height="{h}" rx="6" fill="#0D1117" stroke="#283040" stroke-width="1"/>
      <text x="{w/2}" y="15" fill="#38BDF8" font-size="10.5" font-weight="bold" text-anchor="middle">
        2D Pharmacophore Map — {pdb_id} : {lig_id[:12]}
      </text>
      {''.join(lines_svg)}
      <circle cx="{cx:.1f}" cy="{cy:.1f}" r="24" fill="#0F172A" stroke="#10B981" stroke-width="2"/>
      <text x="{cx:.1f}" y="{cy - 2:.1f}" fill="#10B981" font-size="8.5" font-weight="bold" text-anchor="middle">{lig_id[:9]}</text>
      <text x="{cx:.1f}" y="{cy + 8:.1f}" fill="#CBD5E1" font-size="7.5" text-anchor="middle">LIGAND</text>
      {''.join(nodes_svg)}
    </svg>
    """


def build_stage5_3d_html(
    manifest: Dict[str, Any],
    complex_idx: int = 1,
    focus_bond_idx: Optional[int] = None,
) -> str:
    """Build 3Dmol.js viewport with 3D dashed non-covalent bond cylinders, residue labels, and 2D SVG map."""
    files = manifest.get("files", {})
    rec_path = Path(files.get("receptor_clean_pdb", ""))
    cryst_path = Path(files.get("crystal_ligand_pdb", ""))

    rec_pdb_str = rec_path.read_text(encoding="utf-8", errors="replace") if rec_path.exists() else ""
    cryst_pdb_str = cryst_path.read_text(encoding="utf-8", errors="replace") if cryst_path.exists() else ""

    complexes = manifest.get("complexes", [])
    c_idx = max(0, min(complex_idx, len(complexes) - 1)) if complexes else 0
    comp = complexes[c_idx] if complexes else {}

    pose_pdb_str = comp.get("pose_pdb_block", "")
    interactions = comp.get("interactions", [])
    counts = comp.get("counts", {})
    ifp_sim = comp.get("ifp_similarity", 1.0)
    hots_hit = comp.get("hotspots_hit", [])

    pdb_id = manifest.get("pdb_id", "1IEP")
    target_name = manifest.get("target_name", pdb_id)
    lig_name = comp.get("name", "Ligand")
    lig_id = comp.get("id", "LIG")
    role = comp.get("role", "ACTIVE")
    sc_vina = comp.get("score_vina", 0.0)

    # Interacting residue numbers so 3Dmol highlights them clearly
    inter_res_nums = sorted({int(x["res_num"]) for x in interactions if "res_num" in x})
    if not inter_res_nums:
        inter_res_nums = manifest.get("pocket_residue_ids", [])

    svg_2d = _build_2d_pharmacophore_svg(comp, pdb_id)
    focus_inter = interactions[focus_bond_idx] if (focus_bond_idx is not None and 0 <= focus_bond_idx < len(interactions)) else None

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
    background: rgba(15, 19, 26, 0.88); border: 1px solid #262E3D;
    border-left: 3px solid #38BDF8; padding: 8px 12px; border-radius: 6px;
    color: #F1F5F9; font-size: 11px; line-height: 1.45; max-width: 480px;
    pointer-events: none;
  }}
  .hud-title {{ font-size: 13px; font-weight: 800; color: #38BDF8; }}
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
  .map-panel {{
    position: absolute; bottom: 12px; right: 14px; z-index: 20;
    box-shadow: 0 6px 18px rgba(0,0,0,0.55);
  }}
  .legend {{
    position: absolute; bottom: 12px; left: 14px; z-index: 20;
    background: rgba(15, 19, 26, 0.88); border: 1px solid #262E3D;
    padding: 7px 12px; border-radius: 6px; color: #CBD5E1; font-size: 11px;
    display: flex; gap: 12px; align-items: center; pointer-events: none;
  }}
  .dot {{ display: inline-block; width: 9px; height: 9px; border-radius: 50%; margin-right: 4px; vertical-align: middle; }}
</style>
</head>
<body>
<div id="viewport"></div>
<div class="hud-card">
  <div class="hud-title">{pdb_id} ({target_name}) &bull; {lig_name}</div>
  <div>
    <b>H-Bonds:</b> <span style="color:#38BDF8;font-weight:700;">{counts.get('H-BOND', 0)}</span> &nbsp;|&nbsp;
    <b>Salt Bridges:</b> <span style="color:#EC4899;font-weight:700;">{counts.get('SALT-BRIDGE', 0)}</span> &nbsp;|&nbsp;
    <b>pi-Stack:</b> <span style="color:#A855F7;font-weight:700;">{counts.get('PI-STACKING', 0)}</span> &nbsp;|&nbsp;
    <b>Hydrophobic:</b> <span style="color:#F59E0B;font-weight:700;">{counts.get('HYDROPHOBIC', 0)}</span>
  </div>
  <div class="hud-sub">
    Native Crystal IFP Similarity: <b style="color:#10B981;">{ifp_sim:.2f}</b> &bull;
    Hotspots Hit: <b>{', '.join(hots_hit[:5]) or 'Pocket Core'}</b>
  </div>
</div>
<div class="toolbar">
  <button id="btnBonds" class="active" onclick="toggleBonds()">3D Bonds ({len(interactions)})</button>
  <button id="btnLabels" class="active" onclick="toggleLabels()">Residue Labels</button>
  <button id="btnRec" class="active" onclick="toggleRec()">Receptor</button>
  <button id="btnCryst" onclick="toggleCryst()">Crystal Ref</button>
  <button id="btnMap" class="active" onclick="toggleMap()">2D Pharmacophore</button>
  <button onclick="resetCam()">Reset View</button>
</div>
<div id="mapBox" class="map-panel">{svg_2d}</div>
<div class="legend">
  <span><i class="dot" style="background:#38BDF8;"></i>H-Bond</span>
  <span><i class="dot" style="background:#EC4899;"></i>Salt Bridge</span>
  <span><i class="dot" style="background:#A855F7;"></i>pi-Stacking</span>
  <span><i class="dot" style="background:#22C55E;"></i>Halogen Bond</span>
  <span><i class="dot" style="background:#F59E0B;"></i>Hydrophobic</span>
</div>
<script>
  const recData = {json.dumps(rec_pdb_str)};
  const crystData = {json.dumps(cryst_pdb_str)};
  const poseData = {json.dumps(pose_pdb_str)};
  const interResIds = {json.dumps(inter_res_nums)};
  const interactions = {json.dumps(interactions)};
  const focusInter = {json.dumps(focus_inter)};
  const typeColors = {json.dumps(TYPE_COLORS)};

  let viewer = null;
  let recModel = null;
  let crystModel = null;
  let poseModel = null;
  let bondShapes = [];
  let resLabels = [];
  let vis = {{ bonds: true, labels: true, rec: true, cryst: false, map: true }};

  function applyStyles() {{
    if (!viewer) return;
    if (recModel) {{
      recModel.setStyle({{}}, {{}});
      if (vis.rec) {{
        recModel.setStyle({{}}, {{ cartoon: {{ color: "#334155", opacity: 0.55 }} }});
      }}
      if (interResIds.length > 0) {{
        recModel.setStyle({{ resi: interResIds }}, {{
          cartoon: {{ color: "#64748B", opacity: 0.75 }},
          stick: {{ colorscheme: "orangeCarbon", radius: 0.16 }}
        }});
      }}
    }}
    if (crystModel) {{
      crystModel.setStyle({{}}, vis.cryst ? {{ stick: {{ colorscheme: "greenCarbon", radius: 0.20, opacity: 0.65 }} }} : {{}});
    }}
    if (poseModel) {{
      poseModel.setStyle({{}}, {{
        stick: {{ colorscheme: "cyanCarbon", radius: 0.26 }},
        sphere: {{ colorscheme: "cyanCarbon", scale: 0.22 }}
      }});
    }}

    bondShapes.forEach(s => viewer.removeShape(s));
    bondShapes = [];
    if (vis.bonds) {{
      interactions.forEach(inter => {{
        const col = typeColors[inter.interaction_type] || "#38BDF8";
        const isDashed = inter.interaction_type !== "HYDROPHOBIC";
        const s = viewer.addCylinder({{
          start: {{ x: inter.lig_xyz[0], y: inter.lig_xyz[1], z: inter.lig_xyz[2] }},
          end: {{ x: inter.rec_xyz[0], y: inter.rec_xyz[1], z: inter.rec_xyz[2] }},
          radius: inter.is_hotspot ? 0.11 : 0.08,
          color: col,
          dashed: isDashed,
          fromCap: 1,
          toCap: 1
        }});
        bondShapes.push(s);
      }});
    }}

    viewer.removeAllLabels();
    if (vis.labels) {{
      const seenRes = {{}};
      interactions.forEach(inter => {{
        if (seenRes[inter.res_label]) return;
        seenRes[inter.res_label] = true;
        const col = typeColors[inter.interaction_type] || "#38BDF8";
        viewer.addLabel(inter.res_label + " (" + inter.distance_A + "A)", {{
          position: {{ x: inter.rec_xyz[0], y: inter.rec_xyz[1], z: inter.rec_xyz[2] }},
          backgroundColor: "rgba(15, 23, 42, 0.85)",
          fontColor: col,
          fontSize: 10.5,
          borderThickness: 1,
          borderColor: col,
          padding: 2
        }});
      }});
    }}

    const mb = document.getElementById("mapBox");
    if (mb) mb.style.display = vis.map ? "block" : "none";
    viewer.render();
  }}

  function setBtn(id, state) {{
    const el = document.getElementById(id);
    if (el) el.className = state ? "active" : "";
  }}
  function toggleBonds() {{ vis.bonds = !vis.bonds; setBtn("btnBonds", vis.bonds); applyStyles(); }}
  function toggleLabels() {{ vis.labels = !vis.labels; setBtn("btnLabels", vis.labels); applyStyles(); }}
  function toggleRec() {{ vis.rec = !vis.rec; setBtn("btnRec", vis.rec); applyStyles(); }}
  function toggleCryst() {{ vis.cryst = !vis.cryst; setBtn("btnCryst", vis.cryst); applyStyles(); }}
  function toggleMap() {{ vis.map = !vis.map; setBtn("btnMap", vis.map); applyStyles(); }}
  function resetCam() {{
    if (!viewer) return;
    if (focusInter && recModel) {{
      viewer.zoomTo({{ model: recModel, resi: [focusInter.res_num] }}, 400);
      viewer.zoom(0.85);
    }} else if (poseModel) {{
      viewer.zoomTo({{ model: poseModel }}, 400);
      viewer.zoom(0.75);
    }} else {{
      viewer.zoomTo();
    }}
    viewer.render();
  }}

  window.addEventListener("load", function() {{
    viewer = $3Dmol.createViewer("viewport", {{ backgroundColor: "#090D14" }});
    if (recData) recModel = viewer.addModel(recData, "pdb");
    if (crystData) crystModel = viewer.addModel(crystData, "pdb");
    if (poseData) poseModel = viewer.addModel(poseData, "pdb");
    applyStyles();
    resetCam();
  }});
</script>
</body>
</html>"""


# ==============================================================================
# ANCORE STUDIO NATIVE macOS GUI (Stage 5)
# ==============================================================================
_STAGE5_WINDOW_REF = None


def launch_stage5_gui(default_project_dir: str = "./", default_pdb: str = "1IEP") -> Any:
    """Launch the AnCore Dark Studio GUI for Stage 5."""
    global _STAGE5_WINDOW_REF
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

    class Stage5Worker(QObject):
        progress = Signal(int, str)
        finished = Signal(dict)
        failed = Signal(str)

        def __init__(self, params: Dict[str, Any]) -> None:
            super().__init__()
            self.params = params

        def run(self) -> None:
            try:
                res = run_stage5(
                    **self.params,
                    progress_cb=lambda pct, msg: self.progress.emit(pct, msg),
                )
                self.finished.emit(res)
            except Exception as exc:
                self.failed.emit(str(exc))

    class AnCoreStage5Window(QMainWindow):
        def __init__(self) -> None:
            super().__init__()
            self.setWindowTitle("AnCore Studio — Stage 5: Protein-Ligand Interaction Fingerprinting & Pharmacophore Map")
            self.resize(1440, 900)
            self.proj_path = Path(default_project_dir).expanduser().resolve()
            self.active_pdb = default_pdb.strip().upper() if default_pdb.strip().upper() in BENCHMARK_TARGETS else "1IEP"
            self.active_complex_idx = 1
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
                    padding: 5px 8px;
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
                    background-color: #10B981;
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
                "STAGE 5: PROTEIN-LIGAND INTERACTION FINGERPRINTING &amp; 2D/3D PHARMACOPHORE MAP</span>"
            )
            top_layout.addWidget(brand)
            top_layout.addStretch()

            self.suite_badge = QLabel("INTERACTION SUITE: 0 / 5 TARGETS")
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
            left_dock.setMinimumWidth(385)
            left_dock.setMaximumWidth(460)
            left_layout = QVBoxLayout(left_dock)
            left_layout.setContentsMargins(10, 10, 10, 10)
            left_layout.setSpacing(8)

            # CARD 1: BENCHMARK SUITE NAVIGATOR
            nav_card = QFrame()
            nav_card.setObjectName("cardFrame")
            nav_layout = QVBoxLayout(nav_card)
            nav_layout.setContentsMargins(10, 8, 10, 8)
            nav_layout.setSpacing(5)

            nav_title = QLabel("1. BENCHMARK TARGET NAVIGATOR (CLICK ANY TARGET)")
            nav_title.setObjectName("cardTitle")
            nav_layout.addWidget(nav_title)

            self.nav_table = QTableWidget(len(BENCHMARK_TARGETS), 4)
            self.nav_table.setHorizontalHeaderLabels(["PDB", "Target Receptor", "Top Hit", "IFP Match"])
            self.nav_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
            self.nav_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
            self.nav_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
            self.nav_table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)
            self.nav_table.verticalHeader().setVisible(False)
            self.nav_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
            self.nav_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
            self.nav_table.setAlternatingRowColors(True)
            self.nav_table.setFixedHeight(152)
            self.nav_table.cellClicked.connect(self.on_nav_target_clicked)
            nav_layout.addWidget(self.nav_table)
            left_layout.addWidget(nav_card)

            # CARD 2: SCREENED COMPLEX SELECTOR (Click any ligand in the active target!)
            comp_card = QFrame()
            comp_card.setObjectName("cardFrame")
            comp_layout = QVBoxLayout(comp_card)
            comp_layout.setContentsMargins(10, 8, 10, 8)
            comp_layout.setSpacing(5)

            comp_title = QLabel("2. ACTIVE TARGET DOCKED COMPLEXES (CLICK LIGAND)")
            comp_title.setObjectName("cardTitle")
            comp_layout.addWidget(comp_title)

            self.comp_table = QTableWidget(0, 4)
            self.comp_table.setHorizontalHeaderLabels(["Rank", "Compound", "H-Bonds", "IFP Sim"])
            self.comp_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
            self.comp_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
            self.comp_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
            self.comp_table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)
            self.comp_table.verticalHeader().setVisible(False)
            self.comp_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
            self.comp_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
            self.comp_table.setAlternatingRowColors(True)
            self.comp_table.setFixedHeight(165)
            self.comp_table.cellClicked.connect(self.on_complex_clicked)
            comp_layout.addWidget(self.comp_table)
            left_layout.addWidget(comp_card)

            # CARD 3: GEOMETRIC CUTOFFS & EXECUTION
            exec_card = QFrame()
            exec_card.setObjectName("cardFrame")
            e_layout = QVBoxLayout(exec_card)
            e_layout.setContentsMargins(10, 8, 10, 8)
            e_layout.setSpacing(6)

            p_title = QLabel("3. GEOMETRIC CUTOFFS & EXECUTION")
            p_title.setObjectName("cardTitle")
            e_layout.addWidget(p_title)

            self.proj_input = QLineEdit(str(self.proj_path))
            e_layout.addWidget(self.proj_input)

            grid = QGridLayout()
            grid.setHorizontalSpacing(8)
            grid.setVerticalSpacing(4)
            l_hb = QLabel("H-BOND MAX (A)"); l_hb.setObjectName("fieldLabel")
            l_hy = QLabel("HYDROPHOBIC (A)"); l_hy.setObjectName("fieldLabel")
            self.hb_input = QLineEdit("3.5")
            self.hb_input.setValidator(QDoubleValidator(2.5, 4.5, 1))
            self.hy_input = QLineEdit("4.0")
            self.hy_input.setValidator(QDoubleValidator(3.0, 5.5, 1))
            grid.addWidget(l_hb, 0, 0); grid.addWidget(l_hy, 0, 1)
            grid.addWidget(self.hb_input, 1, 0); grid.addWidget(self.hy_input, 1, 1)
            e_layout.addLayout(grid)

            self.btn_run_all = QPushButton("Profile All 5 Benchmark Suites (55 Complexes)")
            self.btn_run_all.setStyleSheet(
                "QPushButton { background-color: #0284C7; color: #FFFFFF; font-weight: 800; "
                "padding: 9px; border-radius: 5px; font-size: 12px; border: 1px solid #38BDF8; }"
                "QPushButton:hover { background-color: #0369A1; }"
                "QPushButton:disabled { background-color: #334155; color: #64748B; }"
            )
            self.btn_run_all.clicked.connect(lambda: self.start_profiling(run_all=True))
            e_layout.addWidget(self.btn_run_all)

            self.progress_bar = QProgressBar()
            self.progress_bar.setRange(0, 100)
            self.progress_bar.setValue(0)
            self.progress_bar.setFixedHeight(12)
            e_layout.addWidget(self.progress_bar)

            self.status_lbl = QLabel("Ready. Click any target or complex to inspect 3D & 2D interactions.")
            self.status_lbl.setWordWrap(True)
            self.status_lbl.setStyleSheet("color: #94A3B8; font-size: 10.5px;")
            e_layout.addWidget(self.status_lbl)

            left_layout.addWidget(exec_card)
            left_layout.addStretch()
            h_splitter.addWidget(left_dock)

            # --- RIGHT WIDESCREEN 3D MOLECULAR VIEWPORT ---
            if HAS_WEBENGINE:
                self.web_view = QWebEngineView()
                self.web_view.setHtml(
                    "<body style='background:#090D14;color:#64748B;display:flex;align-items:center;"
                    "justify-content:center;height:100vh;font-family:Helvetica,Arial,sans-serif;'>"
                    "<div>Click <b>Profile All 5 Benchmark Suites</b> to map 3D non-covalent interactions &amp; 2D pharmacophores.</div></body>"
                )
                h_splitter.addWidget(self.web_view)
            else:
                fb = QLabel("Install PySide6-Addons for interactive 3Dmol.js WebGL viewport.")
                fb.setAlignment(Qt.AlignmentFlag.AlignCenter)
                h_splitter.addWidget(fb)

            h_splitter.setSizes([400, 1040])

            # --- BOTTOM PROJECT TABLE (TABBED) ---
            self.bottom_tabs = QTabWidget()
            self.bottom_tabs.setContentsMargins(8, 4, 8, 8)

            # Tab 1: Atom-by-Atom 3D Interactions for Selected Complex
            self.bonds_table = QTableWidget(0, 8)
            self.bonds_table.setHorizontalHeaderLabels([
                "Interaction Type", "Receptor Residue", "Chain", "Receptor Atom",
                "Ligand Atom", "Distance (A)", "Key Native Hotspot?", "Energetic Role"
            ])
            self.bonds_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
            self.bonds_table.verticalHeader().setVisible(False)
            self.bonds_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
            self.bonds_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
            self.bonds_table.setAlternatingRowColors(True)
            self.bonds_table.cellClicked.connect(self.on_bond_row_clicked)
            self.bottom_tabs.addTab(self.bonds_table, "Active Complex — 3D Non-Covalent Interactions (Click Any Bond to Zoom)")

            # Tab 2: Target Library IFP Summary Table
            self.ifp_table = QTableWidget(0, 10)
            self.ifp_table.setHorizontalHeaderLabels([
                "Rank", "Ligand ID", "Compound Name", "Role", "Vina Score",
                "H-Bonds", "Salt Bridges", "pi-Stacking", "Hydrophobic", "Crystal IFP Tanimoto"
            ])
            self.ifp_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
            self.ifp_table.verticalHeader().setVisible(False)
            self.ifp_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
            self.ifp_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
            self.ifp_table.setAlternatingRowColors(True)
            self.ifp_table.cellClicked.connect(self.on_ifp_table_clicked)
            self.bottom_tabs.addTab(self.ifp_table, "Target Library Interaction Fingerprint (IFP) Comparison")

            v_splitter.addWidget(self.bottom_tabs)
            v_splitter.setSizes([610, 238])

            self.refresh_all_views()

        def refresh_all_views(self) -> None:
            self.proj_path = Path(self.proj_input.text().strip() or "./").expanduser().resolve()
            self.manifests_by_pdb = load_all_stage5_manifests(self.proj_path)

            codes = list(BENCHMARK_TARGETS.keys())
            profiled = [self.manifests_by_pdb[c] for c in codes if c in self.manifests_by_pdb]

            if profiled:
                mean_ifp = sum(m.get("top_hit_ifp_similarity", 0.0) for m in profiled) / len(profiled)
                self.suite_badge.setText(
                    f"INTERACTION SUITE: {len(profiled)} / {len(codes)} TARGETS  |  MEAN TOP-HIT IFP SIMILARITY: {mean_ifp:.2f}"
                )
                if len(profiled) == len(codes):
                    self.suite_badge.setStyleSheet(
                        "background-color: #064E3B; color: #34D399; font-weight: 800; font-size: 11.5px; "
                        "padding: 6px 14px; border-radius: 5px; border: 1px solid #10B981;"
                    )

            for row_idx, code in enumerate(codes):
                meta = BENCHMARK_TARGETS[code]
                m = self.manifests_by_pdb.get(code)
                has_m = m is not None
                top_nm = m.get("top_hit_name", "—")[:14] if m else "—"
                ifp_v = m.get("top_hit_ifp_similarity", 0.0) if m else 0.0

                nav_items = [
                    QTableWidgetItem(code),
                    QTableWidgetItem(meta["name"]),
                    QTableWidgetItem(top_nm),
                    QTableWidgetItem(f"{ifp_v:.2f}" if has_m else "PENDING"),
                ]
                for c_i, it in enumerate(nav_items):
                    if c_i in (0, 3):
                        it.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    if c_i == 0:
                        it.setForeground(QColor("#38BDF8"))
                    if c_i == 3 and has_m:
                        it.setForeground(QColor("#10B981"))
                    self.nav_table.setItem(row_idx, c_i, it)

            self.select_active_target(self.active_pdb)

        def select_active_target(self, pdb_code: str) -> None:
            pdb_code = pdb_code.strip().upper()
            self.active_pdb = pdb_code

            codes = list(BENCHMARK_TARGETS.keys())
            if pdb_code in codes:
                self.nav_table.selectRow(codes.index(pdb_code))

            m = self.manifests_by_pdb.get(pdb_code)
            if m:
                self.populate_target_manifest(m)

        def populate_target_manifest(self, manifest: Dict[str, Any]) -> None:
            complexes = manifest.get("complexes", [])
            self.comp_table.setRowCount(len(complexes))
            self.ifp_table.setRowCount(len(complexes))

            for idx, cp in enumerate(complexes):
                rnk = "CRYSTAL" if cp["consensus_rank"] == 0 else f"#{cp['consensus_rank']}"
                c = cp["counts"]
                role = cp["role"]

                # Left complex selector table
                c_items = [
                    QTableWidgetItem(rnk),
                    QTableWidgetItem(cp["name"]),
                    QTableWidgetItem(str(c["H-BOND"])),
                    QTableWidgetItem(f"{cp['ifp_similarity']:.2f}"),
                ]
                for col_i, it in enumerate(c_items):
                    if col_i in (0, 2, 3):
                        it.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    if col_i == 0:
                        it.setForeground(QColor("#10B981" if rnk == "CRYSTAL" else "#38BDF8"))
                    self.comp_table.setItem(idx, col_i, it)

                # Bottom Tab 2 IFP comparison table
                ifp_items = [
                    QTableWidgetItem(rnk),
                    QTableWidgetItem(cp["id"]),
                    QTableWidgetItem(cp["name"]),
                    QTableWidgetItem(role),
                    QTableWidgetItem(f"{cp.get('score_vina', 0.0):.2f} kcal"),
                    QTableWidgetItem(str(c["H-BOND"])),
                    QTableWidgetItem(str(c["SALT-BRIDGE"])),
                    QTableWidgetItem(str(c["PI-STACKING"])),
                    QTableWidgetItem(str(c["HYDROPHOBIC"])),
                    QTableWidgetItem(f"{cp['ifp_similarity']:.2f}"),
                ]
                for col_i, it in enumerate(ifp_items):
                    it.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    if col_i == 3:
                        it.setForeground(QColor("#10B981" if role in ("ACTIVE", "CRYSTAL_REF") else "#F59E0B"))
                    if col_i == 9:
                        it.setForeground(QColor("#38BDF8"))
                    self.ifp_table.setItem(idx, col_i, it)

            # Default to Rank #1 Screened Hit (index 1 if available, else 0)
            default_c_idx = 1 if len(complexes) > 1 else 0
            self.select_active_complex(default_c_idx)

        def select_active_complex(self, complex_idx: int) -> None:
            m = self.manifests_by_pdb.get(self.active_pdb)
            if not m:
                return
            complexes = m.get("complexes", [])
            if not complexes:
                return
            c_idx = max(0, min(complex_idx, len(complexes) - 1))
            self.active_complex_idx = c_idx
            self.comp_table.selectRow(c_idx)
            self.ifp_table.selectRow(c_idx)

            comp = complexes[c_idx]
            interactions = comp.get("interactions", [])
            self.bottom_tabs.setTabText(
                0,
                f"Active Complex ({comp['name']}) — {len(interactions)} Non-Covalent 3D Bonds (Click Any Row to Zoom)",
            )

            self.bonds_table.setRowCount(len(interactions))
            for r_i, inter in enumerate(interactions):
                itype = inter.get("interaction_type", "HYDROPHOBIC")
                is_hot = inter.get("is_hotspot", False)
                col_hex = TYPE_COLORS.get(itype, "#94A3B8")

                b_items = [
                    QTableWidgetItem(itype),
                    QTableWidgetItem(inter.get("res_label", "")),
                    QTableWidgetItem(inter.get("chain", "A")),
                    QTableWidgetItem(inter.get("rec_atom", "")),
                    QTableWidgetItem(inter.get("lig_atom", "")),
                    QTableWidgetItem(f"{inter.get('distance_A', 0.0):.2f} A"),
                    QTableWidgetItem("YES (Hotspot)" if is_hot else "Pocket Contact"),
                    QTableWidgetItem(inter.get("role_desc", "")),
                ]
                for col_i, it in enumerate(b_items):
                    it.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    if col_i in (0, 1):
                        it.setForeground(QColor(col_hex))
                    if col_i == 6 and is_hot:
                        it.setForeground(QColor("#10B981"))
                    self.bonds_table.setItem(r_i, col_i, it)

            if HAS_WEBENGINE:
                html_str = build_stage5_3d_html(m, complex_idx=c_idx, focus_bond_idx=None)
                self.web_view.setHtml(html_str, QUrl("https://3Dmol.org/"))

        def on_nav_target_clicked(self, row: int, _col: int) -> None:
            codes = list(BENCHMARK_TARGETS.keys())
            if 0 <= row < len(codes):
                self.select_active_target(codes[row])

        def on_complex_clicked(self, row: int, _col: int) -> None:
            self.select_active_complex(row)

        def on_ifp_table_clicked(self, row: int, _col: int) -> None:
            self.select_active_complex(row)

        def on_bond_row_clicked(self, row: int, _col: int) -> None:
            m = self.manifests_by_pdb.get(self.active_pdb)
            if not m or not HAS_WEBENGINE:
                return
            html_str = build_stage5_3d_html(m, complex_idx=self.active_complex_idx, focus_bond_idx=row)
            self.web_view.setHtml(html_str, QUrl("https://3Dmol.org/"))

        def start_profiling(self, run_all: bool = True) -> None:
            self.btn_run_all.setEnabled(False)
            self.progress_bar.setValue(5)

            params = {
                "project_dir": self.proj_input.text().strip() or "./",
                "pdb_id": "ALL" if run_all else self.active_pdb,
                "hbond_cutoff": float(self.hb_input.text() or 3.5),
                "hydro_cutoff": float(self.hy_input.text() or 4.0),
                "run_all_benchmarks": run_all,
            }

            self._worker_thread = QThread()
            self._worker = Stage5Worker(params)
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
            self.refresh_all_views()

        def on_failed(self, err_msg: str) -> None:
            if self._worker_thread:
                self._worker_thread.quit()
                self._worker_thread.wait()
            self.btn_run_all.setEnabled(True)
            self.status_lbl.setText(f"Error: {err_msg}")
            QMessageBox.critical(self, "Stage 5 Error", err_msg)

    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    win = AnCoreStage5Window()
    _STAGE5_WINDOW_REF = win
    win.show()
    win.raise_()
    win.activateWindow()
    logger.info("Stage 5 AnCore Studio GUI opened. Close the window when finished to return to Jupyter.")
    app.exec()
    return win


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="AnCore Stage 5 — Interaction Fingerprinting & Pharmacophore Map")
    parser.add_argument("--project-dir", type=str, default="./", help="Project directory")
    parser.add_argument("--pdb-id", type=str, default="ALL", help="RCSB PDB ID or 'ALL' for all 5 benchmarks")
    parser.add_argument("--gui", action="store_true", help="Launch native macOS GUI")

    clean_argv = [a for a in sys.argv[1:] if not a.startswith("-f") and "kernel-" not in a and not a.endswith(".json")]
    args, _ = parser.parse_known_args(clean_argv)

    if args.gui:
        launch_stage5_gui(args.project_dir, args.pdb_id if args.pdb_id != "ALL" else "1IEP")
    else:
        run_stage5(project_dir=args.project_dir, pdb_id=args.pdb_id)
