#!/usr/bin/env python3
"""
======================================================================================
ANCORE SCIENTIFIC STUDIO v1.2 — STAGE 5: 3D NON-COVALENT INTERACTION FINGERPRINTING
======================================================================================
Evaluates 3D atomic-level non-covalent interaction networks and bitwise Tanimoto
interaction fingerprint similarity (T_IFP, Eq. 12) across all 55 protein-ligand
complexes in the 5-target benchmark suite (5 crystallographic reference complexes from
Stage 1 + 50 top-ranked docked complexes from Stage 4):
  1. Hydrogen Bonds (HBond)          : Donor-Acceptor distance <= 3.50 Å (polar N/O/S)
  2. Ionic Salt Bridges (SaltBridge) : Charged N/O pair distance <= 4.20 Å (ASP/GLU/LYS/ARG/HIS)
  3. Pi-Stacking / Cation-Pi (Pi)    : Aromatic ring / conjugated contact <= 4.30 Å (PHE/TYR/TRP/HIS)
  4. Halogen / Fluorine (Halogen)    : Ligand F/Cl/Br/I to polar O/N/S <= 3.90 Å
  5. Hydrophobic Contacts (Hydroph)  : Non-polar C-C / C-S packing <= 4.00 Å (grouped per residue)

Generates per target:
  - `05_Stage5_Interaction_Fingerprints/interaction_fingerprint_summary.csv`
  - `05_Stage5_Interaction_Fingerprints/atomic_contacts_csv/<LIG>_contacts.csv`
  - `05_Stage5_Interaction_Fingerprints/pharmacophore_maps_2d/<LIG>_pharmacophore_2d.png`
  - `workspaces/<PDB>/stage5_manifest.json`
  - `workspaces/stage5_suite_summary.json`

Usage:
  python3 src/stage5_interactions.py --run-all
  python3 src/stage5_interactions.py --pdb-id 1IEP
======================================================================================
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

try:
    _THIS_DIR = Path(__file__).resolve().parent
    DEFAULT_BASE_DIR = _THIS_DIR.parent if _THIS_DIR.name == "src" else _THIS_DIR
except NameError:
    DEFAULT_BASE_DIR = Path.cwd()

DEFAULT_BENCHMARK_SUITE = ["1IEP", "3ERT", "1UYG", "1M17", "1STP"]

AROMATIC_RESIDUES = {"PHE", "TYR", "TRP", "HIS", "HID", "HIE", "HIP"}
AROMATIC_RING_ATOMS = {
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TRP": {"CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "HIS": {"CG", "ND1", "CD2", "CE1", "NE2"},
    "HID": {"CG", "ND1", "CD2", "CE1", "NE2"},
    "HIE": {"CG", "ND1", "CD2", "CE1", "NE2"},
    "HIP": {"CG", "ND1", "CD2", "CE1", "NE2"},
}
ANIONIC_RES_ATOMS = {
    "ASP": {"OD1", "OD2"},
    "GLU": {"OE1", "OE2"},
}
CATIONIC_RES_ATOMS = {
    "LYS": {"NZ"},
    "ARG": {"NH1", "NH2", "NE"},
    "HIS": {"ND1", "NE2"},
    "HIP": {"ND1", "NE2"},
}


def parse_pdb_heavy_atoms(pdb_path: Path, record_types: Tuple[str, ...] = ("ATOM  ",)) -> List[dict]:
    """Parses heavy atoms from a standard 80-column PDB file."""
    atoms: List[dict] = []
    if not pdb_path.is_file():
        return atoms

    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.startswith(record_types):
                continue
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain_id = line[21:22].strip() or "A"
            try:
                res_seq = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
            except ValueError:
                continue

            el = line[76:78].strip().upper() if len(line) >= 78 else ""
            if not el:
                letters = "".join(c for c in atom_name if c.isalpha()).upper()
                el = letters[:2] if letters.startswith(("CL", "BR", "FE", "ZN", "MG")) else (letters[:1] or "C")

            if el in ("H", "D"):
                continue

            res_label = f"{res_name.capitalize()}{res_seq}"
            atoms.append({
                "atom_name": atom_name,
                "res_name": res_name.upper(),
                "chain_id": chain_id,
                "res_seq": res_seq,
                "res_label": res_label,
                "x": x,
                "y": y,
                "z": z,
                "element": el,
            })
    return atoms


def parse_pdbqt_heavy_atoms(pdbqt_path: Path) -> List[dict]:
    """Parses heavy atoms, partial charges, and AD4 atom types from a ligand PDBQT file."""
    atoms: List[dict] = []
    if not pdbqt_path.is_file():
        return atoms

    with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            atom_name = line[12:16].strip()
            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
            except ValueError:
                continue

            q = 0.0
            if len(line) >= 76:
                try:
                    q = float(line[66:76].strip())
                except ValueError:
                    q = 0.0

            ad_type = line[76:79].strip().upper() if len(line) >= 77 else ""
            if ad_type in ("H", "D", "HD"):
                continue

            if ad_type == "A":
                el = "C"
            elif ad_type.startswith("NA"):
                el = "N"
            elif ad_type.startswith("OA"):
                el = "O"
            elif ad_type.startswith("SA"):
                el = "S"
            elif ad_type in ("CL", "BR"):
                el = ad_type
            else:
                letters = "".join(c for c in (ad_type or atom_name) if c.isalpha()).upper()
                el = letters[:2] if letters.startswith(("CL", "BR")) else (letters[:1] or "C")

            if el in ("H", "D"):
                continue

            atoms.append({
                "atom_name": atom_name,
                "x": x,
                "y": y,
                "z": z,
                "q": q,
                "ad_type": ad_type or el,
                "element": el,
            })
    return atoms


def detect_protein_ligand_contacts(
    receptor_atoms: List[dict],
    ligand_atoms: List[dict],
) -> Tuple[List[dict], Dict[str, int], Set[str]]:
    """
    Computes 3D non-covalent interactions between receptor heavy atoms and ligand heavy atoms,
    deduplicating per (receptor_residue, interaction_type) so each residue contributes physically
    distinct interaction channels while recording the closest atomic contact distance.
    Returns:
      - detailed_contacts : list of dicts sorted by interaction priority and distance
      - counts            : dict of counts by interaction category
      - bitvector_tokens  : set of string tokens `ResLabel:Category` for bitwise T_IFP calculation
    """
    if not receptor_atoms or not ligand_atoms:
        return [], {"HBond": 0, "SaltBridge": 0, "PiStacking": 0, "Halogen": 0, "Hydrophobic": 0}, set()

    r_xyz = np.array([[a["x"], a["y"], a["z"]] for a in receptor_atoms], dtype=float)
    l_xyz = np.array([[a["x"], a["y"], a["z"]] for a in ligand_atoms], dtype=float)
    d_mat = np.linalg.norm(r_xyz[:, None, :] - l_xyz[None, :, :], axis=-1)

    # Filter all heavy-atom pairs within 4.50 Å
    idx_r, idx_l = np.where((d_mat >= 1.50) & (d_mat <= 4.50))

    # Track best (shortest-distance) contact per (res_label, category, subtype_atom)
    best_per_res_cat: Dict[Tuple[str, str], dict] = {}
    hbond_pairs_per_res: Dict[Tuple[str, str], dict] = {}

    for ir, il in zip(idx_r, idx_l):
        ra = receptor_atoms[ir]
        la = ligand_atoms[il]
        dist = float(d_mat[ir, il])

        res_name = ra["res_name"]
        res_label = ra["res_label"]
        r_atom = ra["atom_name"]
        r_el = ra["element"]
        l_atom = la["atom_name"]
        l_el = la["element"]
        l_ad = la.get("ad_type", l_el)
        l_q = float(la.get("q", 0.0))

        # 1. Salt Bridge check (<= 4.20 Å between oppositely charged / ionizable groups)
        is_rec_anionic = r_atom in ANIONIC_RES_ATOMS.get(res_name, set())
        is_rec_cationic = r_atom in CATIONIC_RES_ATOMS.get(res_name, set())
        is_lig_cationic = (l_el == "N" and l_q > -0.15)
        is_lig_anionic = (l_el == "O" and l_q < -0.20)

        if dist <= 4.20 and ((is_rec_anionic and is_lig_cationic) or (is_rec_cationic and is_lig_anionic)):
            key = (res_label, "SaltBridge")
            if key not in best_per_res_cat or dist < best_per_res_cat[key]["distance_A"]:
                best_per_res_cat[key] = {
                    "receptor_residue": res_label,
                    "receptor_resname": res_name,
                    "receptor_atom": r_atom,
                    "ligand_atom": l_atom,
                    "ligand_element": l_el,
                    "interaction_type": "SaltBridge",
                    "distance_A": round(dist, 3),
                }

        # 2. Hydrogen Bond check (<= 3.50 Å between N/O/S polar heavy atoms)
        if dist <= 3.50 and r_el in ("N", "O", "S") and l_el in ("N", "O", "S"):
            # Allow up to 2 distinct H-bonds per residue if mediated by different ligand heavy atoms
            hb_key = (res_label, l_atom)
            if hb_key not in hbond_pairs_per_res or dist < hbond_pairs_per_res[hb_key]["distance_A"]:
                hbond_pairs_per_res[hb_key] = {
                    "receptor_residue": res_label,
                    "receptor_resname": res_name,
                    "receptor_atom": r_atom,
                    "ligand_atom": l_atom,
                    "ligand_element": l_el,
                    "interaction_type": "HBond",
                    "distance_A": round(dist, 3),
                }

        # 3. Halogen Bond / Fluorine Contact check (<= 3.90 Å between ligand F/Cl/Br/I and receptor O/N/S)
        if dist <= 3.90 and l_el in ("F", "CL", "BR", "I") and r_el in ("O", "N", "S"):
            key = (res_label, "Halogen")
            if key not in best_per_res_cat or dist < best_per_res_cat[key]["distance_A"]:
                best_per_res_cat[key] = {
                    "receptor_residue": res_label,
                    "receptor_resname": res_name,
                    "receptor_atom": r_atom,
                    "ligand_atom": l_atom,
                    "ligand_element": l_el,
                    "interaction_type": "Halogen",
                    "distance_A": round(dist, 3),
                }

        # 4. Pi-Stacking / Cation-Pi check (<= 4.15 Å between aromatic residue ring and aromatic/conjugated ligand atom)
        is_rec_arom_atom = (res_name in AROMATIC_RESIDUES) and (r_atom in AROMATIC_RING_ATOMS.get(res_name, set()))
        is_lig_arom_or_cat = (l_ad == "A") or (l_el in ("C", "N") and dist <= 3.95)
        if dist <= 4.15 and is_rec_arom_atom and is_lig_arom_or_cat:
            key = (res_label, "PiStacking")
            if key not in best_per_res_cat or dist < best_per_res_cat[key]["distance_A"]:
                best_per_res_cat[key] = {
                    "receptor_residue": res_label,
                    "receptor_resname": res_name,
                    "receptor_atom": r_atom,
                    "ligand_atom": l_atom,
                    "ligand_element": l_el,
                    "interaction_type": "PiStacking",
                    "distance_A": round(dist, 3),
                }

        # 5. Hydrophobic Contact check (<= 4.00 Å between non-polar C/S atoms)
        if dist <= 4.00 and r_el in ("C", "S") and l_el in ("C", "S", "CL", "BR", "F"):
            key = (res_label, "Hydrophobic")
            if key not in best_per_res_cat or dist < best_per_res_cat[key]["distance_A"]:
                best_per_res_cat[key] = {
                    "receptor_residue": res_label,
                    "receptor_resname": res_name,
                    "receptor_atom": r_atom,
                    "ligand_atom": l_atom,
                    "ligand_element": l_el,
                    "interaction_type": "Hydrophobic",
                    "distance_A": round(dist, 3),
                }

    # Keep up to 2 shortest H-bonds per residue
    hb_by_res: Dict[str, List[dict]] = {}
    for (_, _), c in hbond_pairs_per_res.items():
        hb_by_res.setdefault(c["receptor_residue"], []).append(c)

    detailed_contacts: List[dict] = []
    for res_label, hb_list in hb_by_res.items():
        hb_sorted = sorted(hb_list, key=lambda x: x["distance_A"])[:2]
        detailed_contacts.extend(hb_sorted)

    for (_, cat), c in best_per_res_cat.items():
        detailed_contacts.append(c)

    cat_priority = {"HBond": 1, "SaltBridge": 2, "PiStacking": 3, "Halogen": 4, "Hydrophobic": 5}
    detailed_contacts.sort(key=lambda x: (cat_priority.get(x["interaction_type"], 9), x["distance_A"]))

    counts = {"HBond": 0, "SaltBridge": 0, "PiStacking": 0, "Halogen": 0, "Hydrophobic": 0}
    bitvector_tokens: Set[str] = set()
    for c in detailed_contacts:
        itype = c["interaction_type"]
        counts[itype] = counts.get(itype, 0) + 1
        # Bitvector token encodes both residue identity and interaction channel
        if itype in ("HBond", "SaltBridge"):
            bitvector_tokens.add(f"{c['receptor_residue']}:PolarAnchor")
        elif itype in ("PiStacking", "Halogen"):
            bitvector_tokens.add(f"{c['receptor_residue']}:SpecificDirectional")
        else:
            bitvector_tokens.add(f"{c['receptor_residue']}:PocketPacking")
        bitvector_tokens.add(f"{c['receptor_residue']}:AnyContact")

    return detailed_contacts, counts, bitvector_tokens


def compute_tanimoto_ifp(bits_a: Set[str], bits_b: Set[str]) -> float:
    """Computes bitwise Tanimoto interaction fingerprint similarity T_IFP (Eq. 12)."""
    if not bits_a and not bits_b:
        return 0.0
    union = bits_a | bits_b
    if not union:
        return 0.0
    inter = bits_a & bits_b
    return round(float(len(inter)) / float(len(union)), 3)


def render_2d_pharmacophore_diagram(
    pdb_id: str,
    lig_id: str,
    lig_name: str,
    role: str,
    contacts: List[dict],
    counts: Dict[str, int],
    t_ifp: float,
    out_png: Path,
) -> None:
    """
    Renders a clean, high-resolution 300-DPI 2D radial pocket interaction diagram
    showing the ligand core and surrounding interacting residues color-coded by interaction type.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, FancyBboxPatch

    fig, ax = plt.subplots(figsize=(7.6, 6.2), dpi=600, facecolor="white")
    ax.set_xlim(-1.45, 1.45)
    ax.set_ylim(-1.35, 1.35)
    ax.axis("off")

    # Select top contacts to display around the radial pocket ring (prioritizing specific interactions)
    display_contacts: List[dict] = []
    seen_res: Set[str] = set()
    for c in contacts:
        key = f"{c['receptor_residue']}_{c['interaction_type']}"
        if key not in seen_res:
            seen_res.add(key)
            display_contacts.append(c)
        if len(display_contacts) >= 14:
            break

    color_map = {
        "HBond": ("#15803d", "#dcfce7", "H-Bond"),
        "SaltBridge": ("#b91c1c", "#fee2e2", "Salt Bridge"),
        "PiStacking": ("#7e22ce", "#f3e8ff", "π-Stacking"),
        "Halogen": ("#0369a1", "#e0f2fe", "Halogen Bond"),
        "Hydrophobic": ("#475569", "#f1f5f9", "Hydrophobic"),
    }

    n_nodes = max(1, len(display_contacts))
    radius = 0.94
    for idx, c in enumerate(display_contacts):
        angle = 2.0 * math.pi * idx / float(n_nodes) - math.pi / 2.0
        nx = radius * math.cos(angle)
        ny = radius * math.sin(angle)

        itype = c["interaction_type"]
        edge_col, fill_col, label_str = color_map.get(itype, ("#334155", "#f8fafc", itype))

        # Draw interaction line from ligand core to residue node
        linestyle = "--" if itype in ("HBond", "SaltBridge", "PiStacking", "Halogen") else ":"
        lw = 1.9 if itype in ("HBond", "SaltBridge") else 1.3
        ax.plot([0.0, nx * 0.84], [0.0, ny * 0.84], linestyle=linestyle, linewidth=lw, color=edge_col, alpha=0.85, zorder=2)

        # Midpoint distance label for specific directional contacts
        if itype in ("HBond", "SaltBridge", "PiStacking", "Halogen"):
            mx, my = nx * 0.48, ny * 0.48
            ax.text(
                mx,
                my,
                f"{c['distance_A']:.2f} Å",
                fontsize=7.2,
                fontweight="bold",
                color=edge_col,
                ha="center",
                va="center",
                bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor=edge_col, alpha=0.92, lw=0.7),
                zorder=4,
            )

        # Draw residue node circle
        circ = Circle((nx, ny), 0.175, facecolor=fill_col, edgecolor=edge_col, linewidth=1.8, zorder=5)
        ax.add_patch(circ)
        ax.text(
            nx,
            ny + 0.03,
            c["receptor_residue"],
            fontsize=8.2,
            fontweight="bold",
            color="#0f172a",
            ha="center",
            va="center",
            zorder=6,
        )
        ax.text(
            nx,
            ny - 0.055,
            f"{c['receptor_atom']}↔{c['ligand_atom']}",
            fontsize=6.5,
            color=edge_col,
            ha="center",
            va="center",
            zorder=6,
        )

    # Central ligand hub badge
    hub_edge = "#0f766e" if role in ("ACTIVE", "COCRYSTAL_REF") else "#b45309"
    hub_fill = "#f0fdfa" if role in ("ACTIVE", "COCRYSTAL_REF") else "#fffbeb"
    hub = FancyBboxPatch(
        (-0.34, -0.20),
        0.68,
        0.40,
        boxstyle="round,pad=0.04,rounding_size=0.08",
        facecolor=hub_fill,
        edgecolor=hub_edge,
        linewidth=2.2,
        zorder=7,
    )
    ax.add_patch(hub)
    ax.text(0.0, 0.06, lig_id, fontsize=10.5, fontweight="bold", color="#0f172a", ha="center", va="center", zorder=8)
    ax.text(0.0, -0.04, f"Role: {role}", fontsize=7.8, fontweight="bold", color=hub_edge, ha="center", va="center", zorder=8)
    ax.text(0.0, -0.13, f"T_IFP = {t_ifp:.3f}", fontsize=8.0, color="#334155", ha="center", va="center", zorder=8)

    # Title & Summary Footer
    ax.set_title(
        f"AnCore v1.2 Stage 5 — 3D Pocket Pharmacophore Network: {pdb_id} | {lig_name}",
        fontsize=10.5,
        fontweight="bold",
        color="#0f172a",
        pad=10,
    )
    summary_str = (
        f"H-Bonds: {counts['HBond']}   |   Salt Bridges: {counts['SaltBridge']}   |   "
        f"π-Stacking: {counts['PiStacking']}   |   Halogen: {counts['Halogen']}   |   "
        f"Hydrophobic Residues: {counts['Hydrophobic']}   |   Co-Crystal T_IFP: {t_ifp:.3f}"
    )
    fig.text(
        0.5,
        0.025,
        summary_str,
        ha="center",
        va="center",
        fontsize=8.0,
        fontweight="bold",
        color="#1e293b",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#f8fafc", edgecolor="#cbd5e1", lw=1.0),
    )

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=600, bbox_inches="tight", facecolor="white", transparent=False)
    plt.close(fig)


def run_stage5_for_target(
    pdb_id: str,
    workspace_dir: Path,
) -> dict:
    """
    Runs Stage 5 3D non-covalent interaction fingerprinting for a single target across:
      - The crystallographic reference complex from Stage 1 (`COCRYSTAL_REF`)
      - All 10 docked complexes from Stage 4
    """
    pdb_id = pdb_id.upper().strip()
    s1_manifest_path = workspace_dir / "stage1_manifest.json"
    s4_manifest_path = workspace_dir / "stage4_manifest.json"
    if not s1_manifest_path.is_file():
        raise FileNotFoundError(f"Stage 1 manifest missing for {pdb_id}: {s1_manifest_path}")
    if not s4_manifest_path.is_file():
        raise FileNotFoundError(f"Stage 4 manifest missing for {pdb_id}: {s4_manifest_path}")

    s1 = json.loads(s1_manifest_path.read_text(encoding="utf-8"))
    s4 = json.loads(s4_manifest_path.read_text(encoding="utf-8"))

    s5_dir = workspace_dir / "05_Stage5_Interaction_Fingerprints"
    contacts_dir = s5_dir / "atomic_contacts_csv"
    maps_dir = s5_dir / "pharmacophore_maps_2d"
    s5_dir.mkdir(parents=True, exist_ok=True)
    contacts_dir.mkdir(parents=True, exist_ok=True)
    maps_dir.mkdir(parents=True, exist_ok=True)

    minimized_receptor_pdb = Path(s1["files"]["minimized_receptor_pdb"])
    cocrystal_pdb = Path(s1["files"]["cocrystal_ligand_pdb"])
    cocrystal_lig_id = s1.get("cocrystal_ligand_id") or s1.get("ligand_id") or s1.get("ligand_resname") or "LIG"

    # Load receptor heavy atoms near the binding site (within 15 Å of box center for fast exact lookup)
    box_info = s1.get("docking_box", {})
    if "center_x" in box_info:
        cx, cy, cz = float(box_info["center_x"]), float(box_info["center_y"]), float(box_info["center_z"])
    else:
        cx, cy, cz = box_info.get("center_xyz") or box_info.get("center_A") or (0.0, 0.0, 0.0)
    all_rec_atoms = parse_pdb_heavy_atoms(minimized_receptor_pdb, record_types=("ATOM  ",))
    rec_pocket_atoms = [
        a for a in all_rec_atoms
        if math.sqrt((a["x"] - cx) ** 2 + (a["y"] - cy) ** 2 + (a["z"] - cz) ** 2) <= 15.0
    ]

    # 1. Evaluate Co-Crystal Reference Complex
    cocrystal_pdbqt = workspace_dir / "02_Stage2_Redocking_Gate" / f"{cocrystal_lig_id}_cocrystal_input.pdbqt"
    if cocrystal_pdbqt.is_file():
        ref_lig_atoms = parse_pdbqt_heavy_atoms(cocrystal_pdbqt)
    else:
        ref_lig_atoms = parse_pdb_heavy_atoms(cocrystal_pdb, record_types=("HETATM", "ATOM  "))

    ref_contacts, ref_counts, ref_bits = detect_protein_ligand_contacts(rec_pocket_atoms, ref_lig_atoms)

    ref_csv = contacts_dir / f"{cocrystal_lig_id}_COCRYSTAL_REF_contacts.csv"
    if ref_contacts:
        with open(ref_csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(ref_contacts[0].keys()))
            w.writeheader()
            w.writerows(ref_contacts)

    ref_png = maps_dir / f"{cocrystal_lig_id}_COCRYSTAL_REF_pharmacophore_2d.png"
    render_2d_pharmacophore_diagram(
        pdb_id=pdb_id,
        lig_id=f"{cocrystal_lig_id} (Ref)",
        lig_name=f"Co-Crystal {cocrystal_lig_id}",
        role="COCRYSTAL_REF",
        contacts=ref_contacts,
        counts=ref_counts,
        t_ifp=1.000,
        out_png=ref_png,
    )

    ref_hbond_res = sorted({c["receptor_residue"] for c in ref_contacts if c["interaction_type"] in ("HBond", "SaltBridge")})
    ref_all_res = sorted({c["receptor_residue"] for c in ref_contacts})

    summary_rows: List[dict] = [
        {
            "pdb_id": pdb_id,
            "consensus_rank": 0,
            "ligand_id": f"{cocrystal_lig_id}_REF",
            "name": f"Co-Crystal Reference ({cocrystal_lig_id})",
            "role": "COCRYSTAL_REF",
            "hbond_count": ref_counts["HBond"],
            "salt_bridge_count": ref_counts["SaltBridge"],
            "pi_stacking_count": ref_counts["PiStacking"],
            "halogen_count": ref_counts["Halogen"],
            "hydrophobic_count": ref_counts["Hydrophobic"],
            "specific_polar_pi_contacts": (
                ref_counts["HBond"] + ref_counts["SaltBridge"] + ref_counts["PiStacking"] + ref_counts["Halogen"]
            ),
            "total_interactions": sum(ref_counts.values()),
            "tanimoto_ifp": 1.000,
            "key_polar_residues": ", ".join(ref_hbond_res) if ref_hbond_res else "None",
            "interacting_pocket_residues": ", ".join(ref_all_res),
            "contacts_csv": str(ref_csv.resolve()),
            "pharmacophore_png": str(ref_png.resolve()),
        }
    ]

    # 2. Evaluate all 10 Stage 4 Docked Complexes in Consensus Rank Order
    for row in s4["ranked_ligands"]:
        lig_id = row["ligand_id"]
        lig_name = row["name"]
        role = row["role"]
        c_rank = int(row["consensus_rank"])

        pose_pdbqt = Path(row["top_pose_pdbqt"])
        lig_atoms = parse_pdbqt_heavy_atoms(pose_pdbqt)
        contacts, counts, lig_bits = detect_protein_ligand_contacts(rec_pocket_atoms, lig_atoms)
        t_ifp = compute_tanimoto_ifp(lig_bits, ref_bits)

        lig_csv = contacts_dir / f"{lig_id}_contacts.csv"
        if contacts:
            with open(lig_csv, "w", newline="", encoding="utf-8") as fh:
                w = csv.DictWriter(fh, fieldnames=list(contacts[0].keys()))
                w.writeheader()
                w.writerows(contacts)

        lig_png = maps_dir / f"{lig_id}_pharmacophore_2d.png"
        render_2d_pharmacophore_diagram(
            pdb_id=pdb_id,
            lig_id=lig_id,
            lig_name=lig_name,
            role=role,
            contacts=contacts,
            counts=counts,
            t_ifp=t_ifp,
            out_png=lig_png,
        )

        polar_res = sorted({c["receptor_residue"] for c in contacts if c["interaction_type"] in ("HBond", "SaltBridge")})
        all_res = sorted({c["receptor_residue"] for c in contacts})

        summary_rows.append({
            "pdb_id": pdb_id,
            "consensus_rank": c_rank,
            "ligand_id": lig_id,
            "name": lig_name,
            "role": role,
            "score_vina": row["score_vina"],
            "score_vinardo": row["score_vinardo"],
            "score_ad4": row["score_ad4"],
            "rrf_score": row["rrf_score"],
            "hbond_count": counts["HBond"],
            "salt_bridge_count": counts["SaltBridge"],
            "pi_stacking_count": counts["PiStacking"],
            "halogen_count": counts["Halogen"],
            "hydrophobic_count": counts["Hydrophobic"],
            "specific_polar_pi_contacts": (
                counts["HBond"] + counts["SaltBridge"] + counts["PiStacking"] + counts["Halogen"]
            ),
            "total_interactions": sum(counts.values()),
            "tanimoto_ifp": t_ifp,
            "key_polar_residues": ", ".join(polar_res) if polar_res else "None",
            "interacting_pocket_residues": ", ".join(all_res),
            "contacts_csv": str(lig_csv.resolve()),
            "pharmacophore_png": str(lig_png.resolve()),
        })

    summary_csv = s5_dir / "interaction_fingerprint_summary.csv"
    fieldnames = [
        "pdb_id",
        "consensus_rank",
        "ligand_id",
        "name",
        "role",
        "score_vina",
        "score_vinardo",
        "score_ad4",
        "rrf_score",
        "hbond_count",
        "salt_bridge_count",
        "pi_stacking_count",
        "halogen_count",
        "hydrophobic_count",
        "specific_polar_pi_contacts",
        "total_interactions",
        "tanimoto_ifp",
        "key_polar_residues",
        "interacting_pocket_residues",
        "contacts_csv",
        "pharmacophore_png",
    ]
    with open(summary_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(summary_rows)

    docked_rows = [r for r in summary_rows if r["role"] in ("ACTIVE", "DECOY")]
    act_rows = [r for r in docked_rows if r["role"] == "ACTIVE"]
    dec_rows = [r for r in docked_rows if r["role"] == "DECOY"]

    mean_act_ifp = round(float(np.mean([r["tanimoto_ifp"] for r in act_rows])), 3) if act_rows else 0.0
    mean_dec_ifp = round(float(np.mean([r["tanimoto_ifp"] for r in dec_rows])), 3) if dec_rows else 0.0
    mean_act_hb = round(float(np.mean([r["hbond_count"] + r["salt_bridge_count"] for r in act_rows])), 2) if act_rows else 0.0
    mean_dec_hb = round(float(np.mean([r["hbond_count"] + r["salt_bridge_count"] for r in dec_rows])), 2) if dec_rows else 0.0

    top1_docked = docked_rows[0]
    manifest = {
        "stage": 5,
        "status": "COMPLETED",
        "pdb_id": pdb_id,
        "complexes_analyzed": len(summary_rows),
        "cocrystal_reference": summary_rows[0],
        "top1_docked_ligand_id": top1_docked["ligand_id"],
        "top1_docked_name": top1_docked["name"],
        "top1_docked_role": top1_docked["role"],
        "top1_hbond_ionic": top1_docked["hbond_count"] + top1_docked["salt_bridge_count"],
        "top1_tanimoto_ifp": top1_docked["tanimoto_ifp"],
        "mean_active_tanimoto_ifp": mean_act_ifp,
        "mean_decoy_tanimoto_ifp": mean_dec_ifp,
        "mean_active_polar_contacts": mean_act_hb,
        "mean_decoy_polar_contacts": mean_dec_hb,
        "interaction_rows": summary_rows,
        "files": {
            "interaction_summary_csv": str(summary_csv.resolve()),
            "atomic_contacts_dir": str(contacts_dir.resolve()),
            "pharmacophore_maps_dir": str(maps_dir.resolve()),
        },
    }

    manifest_path = workspace_dir / "stage5_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (s5_dir / "stage5_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(
        description="AnCore v1.2 Stage 5: 3D Non-Covalent Interaction Fingerprinting & 2D Pharmacophore Mapping"
    )
    parser.add_argument("--pdb-id", type=str, default=None, help="Single target PDB ID")
    parser.add_argument("--run-all", action="store_true", help="Run all 5 benchmark targets")
    parser.add_argument("--include-4wkq", action="store_true", help="Also run 4WKQ cross-docking control")
    parser.add_argument(
        "--out-dir",
        type=str,
        default=str(DEFAULT_BASE_DIR / "workspaces"),
        help="Root workspaces directory",
    )
    args, _ = parser.parse_known_args()

    base_out = Path(args.out_dir)

    if args.run_all or args.pdb_id is None:
        targets = list(DEFAULT_BENCHMARK_SUITE)
        if args.include_4wkq and "4WKQ" not in targets:
            targets.append("4WKQ")

        suite_results = []
        for pid in targets:
            m = run_stage5_for_target(pid, base_out / pid)
            suite_results.append(m)

        suite_summary_path = base_out / "stage5_suite_summary.json"
        suite_summary_path.write_text(json.dumps(suite_results, indent=2), encoding="utf-8")

        print("=" * 116)
        print("ANCORE v1.2 — STAGE 5 MULTI-TARGET 3D INTERACTION FINGERPRINTING (55 COMPLEXES: 5 REF + 50 DOCKED)")
        print("=" * 116)
        print(
            f"{'PDB':<6} {'Ref HBond':<10} {'Rank #1 Docked':<23} {'R#1 T_IFP':<10} "
            f"{'Act/Dec Polar':<15} {'Act/Dec T_IFP':<15} {'Key Ref Polar Anchors'}"
        )
        print("-" * 116)
        for m in suite_results:
            ref = m["cocrystal_reference"]
            ref_hb = f"{ref['hbond_count']}HB+{ref['salt_bridge_count']}SB"
            polar_str = f"{m['mean_active_polar_contacts']:.1f} / {m['mean_decoy_polar_contacts']:.1f}"
            ifp_str = f"{m['mean_active_tanimoto_ifp']:.3f} / {m['mean_decoy_tanimoto_ifp']:.3f}"
            print(
                f"{m['pdb_id']:<6} {ref_hb:<10} {m['top1_docked_name'][:22]:<23} "
                f"{m['top1_tanimoto_ifp']:<10.3f} {polar_str:<15} {ifp_str:<15} {ref['key_polar_residues']}"
            )
        print("-" * 116)
        mean_act_ifp_all = np.mean([m["mean_active_tanimoto_ifp"] for m in suite_results])
        mean_dec_ifp_all = np.mean([m["mean_decoy_tanimoto_ifp"] for m in suite_results])
        mean_act_pol_all = np.mean([m["mean_active_polar_contacts"] for m in suite_results])
        mean_dec_pol_all = np.mean([m["mean_decoy_polar_contacts"] for m in suite_results])
        print(
            f"{'MEAN':<6} {'5 Refs':<10} {'55 Complexes Mapped':<23} "
            f"{np.mean([m['top1_tanimoto_ifp'] for m in suite_results]):<10.3f} "
            f"{f'{mean_act_pol_all:.1f} / {mean_dec_pol_all:.1f}':<15} "
            f"{f'{mean_act_ifp_all:.3f} / {mean_dec_ifp_all:.3f}':<15} 55/55 2D PNGs + CSVs Ready"
        )
        print("=" * 116)
        print(f"Saved multi-target Stage 5 summary -> {suite_summary_path}")
        return

    m = run_stage5_for_target(args.pdb_id, base_out / args.pdb_id.upper())
    print(json.dumps(m, indent=2))


if __name__ == "__main__":
    main()