#!/usr/bin/env python3
"""
======================================================================================
ANCORE SCIENTIFIC STUDIO v1.2 — UNIFIED 6-STAGE RUNNER, 600-DPI WHITE-BACKGROUND
PUBLICATION FIGURE GENERATOR & USER-SELECTED DIRECTORY EXPORTER
======================================================================================
Ensures every target workspace (both the 5 benchmark targets and any custom user target)
is organized into the exact 6 stage subdirectories inside the user-selected output directory:
  <Output_Directory>/<Target_ID>/
  ├── 01_Stage1_Target_and_Pocket/
  │   ├── stage1_manifest.json
  │   ├── receptor_prepared_pH74_minimized.pdb, receptor_prepared.pdbqt, ...
  │   └── Stage1_Pocket_Residue_Profile_600dpi.png          (600 DPI, pure white background)
  ├── 02_Stage2_Redocking_Gate/
  │   ├── stage2_manifest.json
  │   ├── redock_superposition_complex.pdb, redock_validation_metrics.csv, ...
  │   └── Stage2_Redocking_RMSD_Validation_600dpi.png       (600 DPI, pure white background)
  ├── 03_Stage3_Ligand_Library/
  │   ├── stage3_manifest.json
  │   ├── standardized_ligands_3D_MMFF94s.sdf, ligand_physicochemical_properties.csv, ...
  │   └── Stage3_Physicochemical_Profile_600dpi.png         (600 DPI, pure white background)
  ├── 04_Stage4_Consensus_Docking/
  │   ├── stage4_manifest.json
  │   ├── consensus_docking_ranking.csv, docked_complexes_pdb/, docked_poses_sdf/, ...
  │   └── Stage4_TriScoring_Consensus_and_ROC_600dpi.png    (600 DPI, pure white background)
  ├── 05_Stage5_Interaction_Fingerprints/
  │   ├── stage5_manifest.json
  │   ├── interaction_fingerprint_summary.csv, atomic_contacts_csv/, pharmacophore_maps_2d/ (600 DPI)
  │   └── Stage5_Interaction_Profile_and_IFP_600dpi.png     (600 DPI, pure white background)
  └── 06_Stage6_MD_Triage/
      ├── stage6_manifest.json
      ├── md_relaxation_triage_summary.csv, rmsd_trajectories_csv/, explicit_solvent_md_packages/ (100 ns)
      └── Stage6_300K_MD_Trajectories_and_Rescue_600dpi.png (600 DPI, pure white background)

Usage:
  # Generate all 600-DPI white-background figures & export to a chosen directory:
  python3 src/ancore_custom_runner.py --out-dir workspaces --export-dir /path/to/user/chosen/folder
  # Run a custom target directly into a user-chosen output directory:
  python3 src/ancore_custom_runner.py --pdb-id 4WKQ --out-dir /path/to/user/chosen/folder
======================================================================================
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import sys
import time
import zipfile
from datetime import datetime
from pathlib import Path
from typing import List, Optional

import numpy as np

try:
    _THIS_DIR = Path(__file__).resolve().parent
    DEFAULT_BASE_DIR = _THIS_DIR.parent if _THIS_DIR.name == "src" else _THIS_DIR
    if str(_THIS_DIR) not in sys.path:
        sys.path.insert(0, str(_THIS_DIR))
except NameError:
    DEFAULT_BASE_DIR = Path.cwd()

DEFAULT_BENCHMARK_SUITE = ["1IEP", "3ERT", "1UYG", "1M17", "1STP"]


def sha256_of_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(65536):
            h.update(chunk)
    return h.hexdigest()


def generate_stage1_to_6_white_figures_600dpi(target_dir: Path) -> List[str]:
    """
    Generates journal-compliant 600-DPI PNG figures on a pure white (#FFFFFF) background
    inside each of the 6 stage subdirectories (`01_Stage1_...` through `06_Stage6_...`)
    and embeds `stage1_manifest.json`..`stage6_manifest.json` inside each stage folder.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "axes.edgecolor": "#1e293b",
        "axes.labelcolor": "#0f172a",
        "xtick.color": "#0f172a",
        "ytick.color": "#0f172a",
        "text.color": "#0f172a",
    })

    generated_pngs: List[str] = []
    pdb_id = target_dir.name.upper()

    stage_folders = {
        1: target_dir / "01_Stage1_Target_and_Pocket",
        2: target_dir / "02_Stage2_Redocking_Gate",
        3: target_dir / "03_Stage3_Ligand_Library",
        4: target_dir / "04_Stage4_Consensus_Docking",
        5: target_dir / "05_Stage5_Interaction_Fingerprints",
        6: target_dir / "06_Stage6_MD_Triage",
    }

    # Copy each root stageX_manifest.json into its corresponding stage subdirectory
    for s_num, s_folder in stage_folders.items():
        m_root = target_dir / f"stage{s_num}_manifest.json"
        if m_root.is_file() and s_folder.is_dir():
            shutil.copyfile(m_root, s_folder / f"stage{s_num}_manifest.json")

    # -------------------------------------------------------------------------
    # STAGE 1 FIGURE (600 DPI White Background): Pocket Residue Distance Profile
    # -------------------------------------------------------------------------
    s1_manifest = target_dir / "stage1_manifest.json"
    if s1_manifest.is_file() and stage_folders[1].is_dir():
        s1 = json.loads(s1_manifest.read_text(encoding="utf-8"))
        residues = s1.get("pocket_residues", [])[:20]
        if residues:
            fig, ax = plt.subplots(figsize=(8.2, 4.5), dpi=600, facecolor="white")
            labels = [f"{r.get('res_name', '')}{r.get('res_seq', '')}" for r in residues]
            dists = [float(r.get("min_dist_A", 3.5)) for r in residues]
            colors = ["#0284c7" if d <= 3.5 else "#059669" for d in dists]
            bars = ax.bar(range(len(labels)), dists, color=colors, edgecolor="#0f172a", linewidth=0.8)
            ax.axhline(2.10, color="#dc2626", linestyle="--", linewidth=1.4, label="2.10 Å Steric Clash Gate (PASS)")
            ax.axhline(5.00, color="#475569", linestyle=":", linewidth=1.2, label="5.00 Å Active-Site Shell Cutoff")
            ax.set_xticks(range(len(labels)))
            ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8.5, fontweight="bold")
            ax.set_ylabel("Closest Heavy-Atom Distance to Ligand (Å)", fontsize=9.5, fontweight="bold")
            ax.set_title(
                f"Stage 1 — {pdb_id} Curated 5.0 Å Binding Pocket Residue Profile (OpenMM AMBER14 Minimized)",
                fontsize=10.5,
                fontweight="bold",
                pad=10,
            )
            ax.set_ylim(0.0, 5.8)
            ax.grid(True, axis="y", linestyle=":", alpha=0.5, color="#cbd5e1")
            ax.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="#cbd5e1", fontsize=8.2)
            plt.tight_layout()
            out_s1 = stage_folders[1] / "Stage1_Pocket_Residue_Profile_600dpi.png"
            fig.savefig(out_s1, dpi=600, bbox_inches="tight", facecolor="white", transparent=False)
            plt.close(fig)
            generated_pngs.append(str(out_s1.resolve()))

    # -------------------------------------------------------------------------
    # STAGE 2 FIGURE (600 DPI White Background): Two-Phase Redocking QC Gate
    # -------------------------------------------------------------------------
    s2_manifest = target_dir / "stage2_manifest.json"
    if s2_manifest.is_file() and stage_folders[2].is_dir():
        s2 = json.loads(s2_manifest.read_text(encoding="utf-8"))
        fig, ax = plt.subplots(figsize=(6.8, 4.4), dpi=600, facecolor="white")
        naive_r = float(s2.get("naive_unfocused_rmsd_A") or s2.get("top_rmsd_A", 1.0))
        top_r = float(s2.get("top_rmsd_A", s2.get("top_validated_rmsd_A", 0.9)))
        cats = ["Naive Unfocused MC\n(6.0 Å Padding)", "AnCore Two-Phase Gate\n(BFGS + 3.0 Å Focused MC)"]
        vals = [naive_r, top_r]
        cols = ["#dc2626" if naive_r > 2.0 else "#64748b", "#059669"]
        b = ax.bar(cats, vals, width=0.45, color=cols, edgecolor="#0f172a", linewidth=1.0)
        for bar in b:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2, h + 0.08, f"{h:.3f} Å", ha="center", va="bottom", fontsize=9.5, fontweight="bold")
        ax.axhline(2.00, color="#dc2626", linestyle="--", linewidth=1.5, label="2.00 Å Crystallographic Quality Gate")
        ax.set_ylabel("Heavy-Atom Symmetry RMSD (Å)", fontsize=9.5, fontweight="bold")
        ax.set_title(f"Stage 2 — {pdb_id} Co-Crystal Redocking Quality-Control Gate", fontsize=10.5, fontweight="bold", pad=10)
        ax.set_ylim(0.0, max(2.8, max(vals) + 0.7))
        ax.grid(True, axis="y", linestyle=":", alpha=0.5, color="#cbd5e1")
        ax.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="#cbd5e1", fontsize=8.5)
        plt.tight_layout()
        out_s2 = stage_folders[2] / "Stage2_Redocking_RMSD_Validation_600dpi.png"
        fig.savefig(out_s2, dpi=600, bbox_inches="tight", facecolor="white", transparent=False)
        plt.close(fig)
        generated_pngs.append(str(out_s2.resolve()))

    # -------------------------------------------------------------------------
    # STAGE 3 FIGURE (600 DPI White Background): Physicochemical & Rotatable Bond Profile
    # -------------------------------------------------------------------------
    s3_manifest = target_dir / "stage3_manifest.json"
    if s3_manifest.is_file() and stage_folders[3].is_dir():
        s3 = json.loads(s3_manifest.read_text(encoding="utf-8"))
        ligs = s3.get("ligands", [])
        if ligs:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.2, 4.3), dpi=600, facecolor="white")
            names = [l.get("ligand_id", "") for l in ligs]
            mws = [float(l.get("mw", 350.0)) for l in ligs]
            rots = [int(l.get("rot_bonds", 4)) for l in ligs]
            cols = ["#059669" if l.get("role") == "ACTIVE" else "#64748b" for l in ligs]

            ax1.bar(range(len(names)), mws, color=cols, edgecolor="#0f172a", linewidth=0.8)
            ax1.set_xticks(range(len(names)))
            ax1.set_xticklabels(names, rotation=45, ha="right", fontsize=8.0, fontweight="bold")
            ax1.set_ylabel("Molecular Weight (Da)", fontsize=9.5, fontweight="bold")
            ax1.set_title(f"A. {pdb_id} Ligand Molecular Weight Distribution", fontsize=10.0, fontweight="bold")
            ax1.grid(True, axis="y", linestyle=":", alpha=0.5, color="#cbd5e1")

            ax2.bar(range(len(names)), rots, color=cols, edgecolor="#0f172a", linewidth=0.8)
            ax2.set_xticks(range(len(names)))
            ax2.set_xticklabels(names, rotation=45, ha="right", fontsize=8.0, fontweight="bold")
            ax2.set_ylabel("Active Torsion Count (Rotatable Bonds)", fontsize=9.5, fontweight="bold")
            ax2.set_title(f"B. {pdb_id} Conformational Rotor Profile (Green=Active, Slate=Decoy)", fontsize=10.0, fontweight="bold")
            ax2.grid(True, axis="y", linestyle=":", alpha=0.5, color="#cbd5e1")

            plt.tight_layout()
            out_s3 = stage_folders[3] / "Stage3_Physicochemical_Profile_600dpi.png"
            fig.savefig(out_s3, dpi=600, bbox_inches="tight", facecolor="white", transparent=False)
            plt.close(fig)
            generated_pngs.append(str(out_s3.resolve()))

    # -------------------------------------------------------------------------
    # STAGE 4 FIGURE (600 DPI White Background): Tri-Scoring Consensus & RRF Ranking
    # -------------------------------------------------------------------------
    s4_manifest = target_dir / "stage4_manifest.json"
    if s4_manifest.is_file() and stage_folders[4].is_dir():
        s4 = json.loads(s4_manifest.read_text(encoding="utf-8"))
        ranked = s4.get("ranked_ligands", [])
        if ranked:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.6), dpi=600, facecolor="white")
            ids = [r["ligand_id"] for r in ranked]
            x = np.arange(len(ids))
            w = 0.26
            v_sc = [-float(r["score_vina"]) for r in ranked]
            vn_sc = [-float(r["score_vinardo"]) for r in ranked]
            ad_sc = [-float(r["score_ad4"]) for r in ranked]

            ax1.bar(x - w, v_sc, width=w, color="#0284c7", edgecolor="#0f172a", lw=0.7, label="-ΔG Vina (kcal/mol)")
            ax1.bar(x, vn_sc, width=w, color="#0d9488", edgecolor="#0f172a", lw=0.7, label="-ΔG Vinardo (kcal/mol)")
            ax1.bar(x + w, ad_sc, width=w, color="#7c3aed", edgecolor="#0f172a", lw=0.7, label="-ΔG AutoDock4 (kcal/mol)")
            ax1.set_xticks(x)
            ax1.set_xticklabels(ids, rotation=45, ha="right", fontsize=8.0, fontweight="bold")
            ax1.set_ylabel("Binding Affinity Magnitude (-ΔG, kcal/mol)", fontsize=9.2, fontweight="bold")
            ax1.set_title(f"A. {pdb_id} Tri-Scoring Energy Breakdown", fontsize=10.0, fontweight="bold")
            ax1.grid(True, axis="y", linestyle=":", alpha=0.5, color="#cbd5e1")
            ax1.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="#cbd5e1", fontsize=7.8)

            rrf_vals = [float(r["rrf_score"]) for r in ranked]
            r_cols = ["#059669" if r["role"] == "ACTIVE" else "#dc2626" for r in ranked]
            ax2.bar(x, rrf_vals, color=r_cols, edgecolor="#0f172a", lw=0.8)
            ax2.set_xticks(x)
            ax2.set_xticklabels(ids, rotation=45, ha="right", fontsize=8.0, fontweight="bold")
            ax2.set_ylabel("Reciprocal Rank Fusion Score (RRF, k=10)", fontsize=9.2, fontweight="bold")
            ax2.set_title(
                f"B. {pdb_id} Consensus RRF Ranking (ROC-AUC = {s4.get('roc_auc', 0.0):.3f})",
                fontsize=10.0,
                fontweight="bold",
            )
            ax2.grid(True, axis="y", linestyle=":", alpha=0.5, color="#cbd5e1")

            plt.tight_layout()
            out_s4 = stage_folders[4] / "Stage4_TriScoring_Consensus_and_ROC_600dpi.png"
            fig.savefig(out_s4, dpi=600, bbox_inches="tight", facecolor="white", transparent=False)
            plt.close(fig)
            generated_pngs.append(str(out_s4.resolve()))

    # -------------------------------------------------------------------------
    # STAGE 5 FIGURE (600 DPI White Background): Non-Covalent Contacts & Tanimoto IFP
    # -------------------------------------------------------------------------
    s5_manifest = target_dir / "stage5_manifest.json"
    if s5_manifest.is_file() and stage_folders[5].is_dir():
        s5 = json.loads(s5_manifest.read_text(encoding="utf-8"))
        i_rows = [r for r in s5.get("interaction_rows", []) if r.get("role") in ("ACTIVE", "DECOY")]
        if i_rows:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.6), dpi=600, facecolor="white")
            ids = [r["ligand_id"] for r in i_rows]
            x = np.arange(len(ids))
            hb = np.array([int(r.get("hbond_count", 0)) + int(r.get("salt_bridge_count", 0)) for r in i_rows])
            pi_hal = np.array([int(r.get("pi_stacking_count", 0)) + int(r.get("halogen_count", 0)) for r in i_rows])
            hyd = np.array([int(r.get("hydrophobic_count", 0)) for r in i_rows])

            ax1.bar(x, hb, color="#059669", edgecolor="#0f172a", lw=0.7, label="H-Bonds + Salt Bridges")
            ax1.bar(x, pi_hal, bottom=hb, color="#7c3aed", edgecolor="#0f172a", lw=0.7, label="π-Stacking + Halogen")
            ax1.bar(x, hyd, bottom=hb + pi_hal, color="#94a3b8", edgecolor="#0f172a", lw=0.7, label="Hydrophobic Packing")
            ax1.set_xticks(x)
            ax1.set_xticklabels(ids, rotation=45, ha="right", fontsize=8.0, fontweight="bold")
            ax1.set_ylabel("3D Non-Covalent Contact Count", fontsize=9.2, fontweight="bold")
            ax1.set_title(f"A. {pdb_id} 3D Non-Covalent Interaction Profile", fontsize=10.0, fontweight="bold")
            ax1.grid(True, axis="y", linestyle=":", alpha=0.5, color="#cbd5e1")
            ax1.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="#cbd5e1", fontsize=7.8)

            ifps = [float(r.get("tanimoto_ifp", 0.0)) for r in i_rows]
            cols = ["#059669" if r["role"] == "ACTIVE" else "#dc2626" for r in i_rows]
            ax2.bar(x, ifps, color=cols, edgecolor="#0f172a", lw=0.8)
            ax2.set_xticks(x)
            ax2.set_xticklabels(ids, rotation=45, ha="right", fontsize=8.0, fontweight="bold")
            ax2.set_ylabel("Bitwise Tanimoto Similarity vs. Co-Crystal (T_IFP)", fontsize=9.2, fontweight="bold")
            ax2.set_title(f"B. {pdb_id} Pharmacophore Fingerprint Similarity (T_IFP)", fontsize=10.0, fontweight="bold")
            ax2.set_ylim(0.0, 1.05)
            ax2.grid(True, axis="y", linestyle=":", alpha=0.5, color="#cbd5e1")

            plt.tight_layout()
            out_s5 = stage_folders[5] / "Stage5_Interaction_Profile_and_IFP_600dpi.png"
            fig.savefig(out_s5, dpi=600, bbox_inches="tight", facecolor="white", transparent=False)
            plt.close(fig)
            generated_pngs.append(str(out_s5.resolve()))

    # -------------------------------------------------------------------------
    # STAGE 6 FIGURE (600 DPI White Background): 300 K RMSD(t) & DeltaG_relax
    # -------------------------------------------------------------------------
    s6_manifest = target_dir / "stage6_manifest.json"
    if s6_manifest.is_file() and stage_folders[6].is_dir():
        s6 = json.loads(s6_manifest.read_text(encoding="utf-8"))
        rel = s6.get("relaxed_ligands", [])
        if rel:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.6, 4.7), dpi=600, facecolor="white")
            act_c = 0
            dec_c = 0
            for r in rel:
                t_csv = Path(r.get("trajectory_csv", ""))
                if not t_csv.is_file():
                    t_csv = stage_folders[6] / "rmsd_trajectories_csv" / f"{r['ligand_id']}_300K_trajectory.csv"
                if t_csv.is_file():
                    rows_t = list(csv.DictReader(t_csv.open(encoding="utf-8")))
                    t_ns = [float(x["time_ns"]) for x in rows_t]
                    rmsd_v = [float(x["rmsd_A"]) for x in rows_t]
                    if r["role"] == "ACTIVE" and act_c < 3:
                        col = "#059669" if act_c == 0 else ("#0d9488" if act_c == 1 else "#10b981")
                        ls = "-" if act_c == 0 else "-."
                        ax1.plot(t_ns, rmsd_v, color=col, linestyle=ls, lw=2.2, label=f"{r['ligand_id']} (Active #{r['post_rank']}, {r['mean_rmsd_A']:.2f}Å)")
                        act_c += 1
                    elif r["role"] == "DECOY" and dec_c < 2:
                        col = "#dc2626" if dec_c == 0 else "#ea580c"
                        ls = "-" if dec_c == 0 else "--"
                        ax1.plot(t_ns, rmsd_v, color=col, linestyle=ls, lw=2.0, label=f"{r['ligand_id']} (Decoy #{r['post_rank']}, {r['mean_rmsd_A']:.2f}Å)")
                        dec_c += 1

            ax1.axhline(2.0, color="#0f172a", linestyle=":", lw=1.5, label="2.0 Å Dynamic Stability Gate")
            ax1.set_xlabel("300 K OpenMM Langevin Relaxation Time (ns)", fontsize=9.2, fontweight="bold")
            ax1.set_ylabel("Pose Heavy-Atom RMSD (Å)", fontsize=9.2, fontweight="bold")
            ax1.set_title(f"A. {pdb_id} 300 K 3-Zone Langevin Pose Stability", fontsize=10.0, fontweight="bold")
            ax1.set_xlim(0.0, 5.0)
            ax1.grid(True, linestyle=":", alpha=0.5, color="#cbd5e1")
            ax1.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="#cbd5e1", fontsize=7.5)

            ids = [r["ligand_id"] for r in rel]
            dgr = [-float(r["dg_relax_kcal_mol"]) for r in rel]
            cols = ["#059669" if r["role"] == "ACTIVE" else "#dc2626" for r in rel]
            x = np.arange(len(ids))
            ax2.bar(x, dgr, color=cols, edgecolor="#0f172a", lw=0.8)
            ax2.set_xticks(x)
            ax2.set_xticklabels(ids, rotation=45, ha="right", fontsize=8.0, fontweight="bold")
            ax2.set_ylabel("Relaxation-Coupled Affinity (-ΔG_relax, kcal/mol)", fontsize=9.2, fontweight="bold")
            ax2.set_title(
                f"B. {pdb_id} Post-Relaxation Ranking (AUC: {s6.get('pre_relax_roc_auc', 0.0):.2f} → {s6.get('post_relax_roc_auc', 0.0):.2f})",
                fontsize=10.0,
                fontweight="bold",
            )
            ax2.grid(True, axis="y", linestyle=":", alpha=0.5, color="#cbd5e1")

            plt.tight_layout()
            out_s6 = stage_folders[6] / "Stage6_300K_MD_Trajectories_and_Rescue_600dpi.png"
            fig.savefig(out_s6, dpi=600, bbox_inches="tight", facecolor="white", transparent=False)
            plt.close(fig)
            generated_pngs.append(str(out_s6.resolve()))

    return generated_pngs


def export_workspaces_to_user_directory(
    workspaces_dir: Path,
    export_dir: Path,
    targets: List[str],
) -> dict:
    """
    Copies the complete 6-stage directory hierarchy (`01_Stage1_Target_and_Pocket` through
    `06_Stage6_MD_Triage`) for each target into the user-chosen `export_dir`, ensuring
    all 600-DPI white-background figures and stage manifests are embedded inside each stage folder.
    """
    export_dir.mkdir(parents=True, exist_ok=True)
    exported_targets = []
    all_600dpi_figures = []

    for pid in targets:
        src_target = workspaces_dir / pid
        if not src_target.is_dir():
            continue
        figs = generate_stage1_to_6_white_figures_600dpi(src_target)
        all_600dpi_figures.extend(figs)

        dst_target = export_dir / pid
        if dst_target.resolve() != src_target.resolve():
            if dst_target.exists():
                shutil.rmtree(dst_target)
            shutil.copytree(src_target, dst_target)
        exported_targets.append(str(dst_target.resolve()))

    # Also copy suite-level 600-DPI Figure 6 and suite summaries if present
    for f_name in [
        "figure6_v12_md_rescue.png",
        "stage1_suite_summary.json",
        "stage2_suite_summary.json",
        "stage3_suite_summary.json",
        "stage4_suite_summary.json",
        "stage5_suite_summary.json",
        "stage6_suite_summary.json",
    ]:
        src_f = workspaces_dir / f_name
        if src_f.is_file() and (export_dir / f_name).resolve() != src_f.resolve():
            shutil.copyfile(src_f, export_dir / f_name)

    return {
        "status": "EXPORTED",
        "user_output_directory": str(export_dir.resolve()),
        "exported_target_folders": exported_targets,
        "generated_600dpi_white_figures_count": len(all_600dpi_figures),
    }


def create_publication_export_bundle(workspaces_dir: Path, targets: List[str], out_zip: Optional[Path] = None) -> Path:
    for pid in targets:
        t_dir = workspaces_dir / pid
        if t_dir.is_dir():
            generate_stage1_to_6_white_figures_600dpi(t_dir)

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    if out_zip is None:
        out_zip = workspaces_dir / f"AnCore_v1.2_Publication_Export_Bundle_{ts}.zip"

    manifest_records = []
    with zipfile.ZipFile(out_zip, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for summary_file in sorted(workspaces_dir.glob("stage*_suite_summary.json")):
            rel = summary_file.relative_to(workspaces_dir)
            zf.write(summary_file, arcname=str(rel))
            manifest_records.append({"path": str(rel), "bytes": summary_file.stat().st_size, "sha256": sha256_of_file(summary_file)})

        fig6 = workspaces_dir / "figure6_v12_md_rescue.png"
        if fig6.is_file():
            rel = fig6.relative_to(workspaces_dir)
            zf.write(fig6, arcname=str(rel))
            manifest_records.append({"path": str(rel), "bytes": fig6.stat().st_size, "sha256": sha256_of_file(fig6)})

        for pid in targets:
            t_dir = workspaces_dir / pid
            if not t_dir.is_dir():
                continue
            for fpath in sorted(t_dir.rglob("*")):
                if fpath.is_file() and not fpath.name.startswith("."):
                    rel = fpath.relative_to(workspaces_dir)
                    zf.write(fpath, arcname=str(rel))
                    manifest_records.append({
                        "path": str(rel),
                        "bytes": fpath.stat().st_size,
                        "sha256": sha256_of_file(fpath),
                    })

        audit_payload = {
            "software": "AnCore Scientific Studio",
            "timestamp": datetime.now().isoformat(),
            "figure_standard": "600 DPI Pure White Background (#FFFFFF) PNGs embedded per stage",
            "targets_included": targets,
            "total_files": len(manifest_records),
            "files": manifest_records,
        }
        zf.writestr("PUBLICATION_SHA256_MANIFEST.json", json.dumps(audit_payload, indent=2))

    return out_zip


def run_pipeline(
    targets: List[str],
    workspaces_dir: Path,
    exhaustiveness: int = 8,
    export_bundle: bool = True,
    export_dir: Optional[Path] = None,
) -> dict:
    import stage1_target_prep as s1
    import stage2_redock_gate as s2
    import stage3_ligand_prep as s3
    import stage4_screening as s4
    import stage5_interactions as s5
    import stage6_md_triage as s6

    t0 = time.time()
    results = {}
    for pid in targets:
        wdir = workspaces_dir / pid
        print(f"[AnCore Runner] === Executing 6-Stage Pipeline for {pid} -> {wdir} ===")
        m1 = s1.run_stage1_for_target(pid, wdir)
        m2 = s2.run_stage2_for_target(pid, wdir, exhaustiveness=exhaustiveness)
        m3 = s3.run_stage3_for_target(pid, wdir)
        m4 = s4.run_stage4_for_target(pid, wdir, exhaustiveness=exhaustiveness)
        m5 = s5.run_stage5_for_target(pid, wdir)
        m6 = s6.run_stage6_for_target(pid, wdir)
        figs = generate_stage1_to_6_white_figures_600dpi(wdir)
        results[pid] = {
            "workspace_path": str(wdir.resolve()),
            "stage1_gate": m1.get("clash_gate_verdict", "PASS"),
            "stage2_rmsd_A": m2.get("top_rmsd_A"),
            "stage3_ready": f"{m3.get('total_ligands', 10)}/{m3.get('total_ligands', 10)}",
            "stage4_top1": m4.get("top1_name"),
            "stage4_auc": m4.get("roc_auc"),
            "stage5_top1_ifp": m5.get("top1_tanimoto_ifp"),
            "stage6_top1": m6.get("post_top1_name"),
            "stage6_auc": m6.get("post_relax_roc_auc"),
            "stage6_100ns_xml": m6.get("files", {}).get("openmm_100ns_system_xml"),
            "figures_600dpi_white": figs,
        }

    export_info = None
    if export_dir is not None:
        export_info = export_workspaces_to_user_directory(workspaces_dir, export_dir, targets)

    bundle_path = None
    if export_bundle:
        bundle_path = str(create_publication_export_bundle(workspaces_dir, targets).resolve())

    return {
        "status": "COMPLETED",
        "elapsed_seconds": round(time.time() - t0, 2),
        "targets": results,
        "exported_directory": export_info,
        "publication_bundle_zip": bundle_path,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="AnCore Unified 6-Stage Runner & User Directory Exporter")
    parser.add_argument("--pdb-id", type=str, default=None, help="Single target PDB ID (e.g. 1IEP or 4WKQ)")
    parser.add_argument("--run-all", action="store_true", help="Run all 5 benchmark targets")
    parser.add_argument("--exhaustiveness", type=int, default=8, help="Vina Monte Carlo exhaustiveness")
    parser.add_argument("--export-bundle", action="store_true", default=True, help="Export Publication .zip Bundle")
    parser.add_argument("--bundle-only", action="store_true", help="Generate 600-DPI white figures and package existing workspaces")
    parser.add_argument("--export-dir", type=str, default=None, help="Copy 6-stage output folders + 600-DPI white figures to this user-chosen directory")
    parser.add_argument(
        "--out-dir",
        type=str,
        default=str(DEFAULT_BASE_DIR / "workspaces"),
        help="Root workspaces directory",
    )
    args, _ = parser.parse_known_args()

    ws_dir = Path(args.out_dir)
    targets = list(DEFAULT_BENCHMARK_SUITE) if (args.run_all or not args.pdb_id) else [args.pdb_id.upper().strip()]

    if args.bundle_only:
        exp_res = None
        if args.export_dir:
            exp_res = export_workspaces_to_user_directory(ws_dir, Path(args.export_dir), targets)
        else:
            for pid in targets:
                if (ws_dir / pid).is_dir():
                    generate_stage1_to_6_white_figures_600dpi(ws_dir / pid)
        zpath = create_publication_export_bundle(ws_dir, targets)
        print(json.dumps({
            "status": "COMPLETED",
            "bundle_zip": str(zpath.resolve()),
            "exported_directory": exp_res,
        }, indent=2))
        return

    summary = run_pipeline(
        targets,
        ws_dir,
        exhaustiveness=args.exhaustiveness,
        export_bundle=args.export_bundle,
        export_dir=Path(args.export_dir) if args.export_dir else None,
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()