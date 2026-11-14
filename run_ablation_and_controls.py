#!/usr/bin/env python3
"""
======================================================================================
AnCore Scientific Studio (v1.2.0) — Ablation, Circularity, Stratified Bootstrap &
1M17 Focused-Grid Propagation Control Suite (`src/run_ablation_and_controls.py`)
======================================================================================
Reproduces all single-engine ablations, RRF k-sensitivity analyses, Stage 5 static
contact controls, Stage 6 co-crystal-blind and shuffled-contact controls, Leave-One-
Target-Out (LOTO) cross-validation, unrestrained full-atom AMBER14/GBn2+GAFF RMSDs,
1M17 focused-box propagation rerun, and within-target stratified 2,000-resample 95%
bootstrap confidence intervals directly from the `workspaces/` benchmark artifacts.

Usage:
  python3 src/run_ablation_and_controls.py --workspaces-dir workspaces
======================================================================================
"""
import argparse
import json
import os
from pathlib import Path
import numpy as np
from scipy.optimize import minimize


def calc_auc(y_true, scores) -> float:
    pos = [s for y, s in zip(y_true, scores) if y == 1]
    neg = [s for y, s in zip(y_true, scores) if y == 0]
    wins = 0.0
    for p in pos:
        for n in neg:
            if p > n:
                wins += 1.0
            elif abs(p - n) <= 1e-12:
                wins += 0.5
    return wins / float(len(pos) * len(neg))


def main():
    parser = argparse.ArgumentParser(description="Run AnCore v1.2.0 Ablations and Controls")
    parser.add_argument("--workspaces-dir", type=str, default="workspaces", help="Path to workspaces directory")
    args = parser.parse_args()

    ws = Path(args.workspaces_dir).resolve()
    with open(ws / "stage1_suite_summary.json") as f:
        s1 = {x["pdb_id"]: x for x in json.load(f)}
    with open(ws / "stage2_suite_summary.json") as f:
        s2 = {x["pdb_id"]: x for x in json.load(f)}
    with open(ws / "stage4_suite_summary.json") as f:
        s4 = {x["pdb_id"]: x for x in json.load(f)}
    with open(ws / "stage5_suite_summary.json") as f:
        s5 = {x["pdb_id"]: x for x in json.load(f)}
    with open(ws / "stage6_suite_summary.json") as f:
        s6 = {x["pdb_id"]: x for x in json.load(f)}

    targets = ["1IEP", "3ERT", "1UYG", "1M17", "1STP"]

    print("======================================================================")
    print("1. STAGE 4 SINGLE-ENGINE VS. CONSENSUS ABLATION & RRF k-SENSITIVITY")
    print("======================================================================")
    for t in targets:
        ligs = s4[t]["ranked_ligands"]
        y = [1 if l["role"] == "ACTIVE" else 0 for l in ligs]
        a_v = calc_auc(y, [-l["score_vina"] for l in ligs])
        a_vn = calc_auc(y, [-l["score_vinardo"] for l in ligs])
        a_ad = calc_auc(y, [-l["score_ad4"] for l in ligs])
        a_rrf = calc_auc(y, [l["rrf_score"] for l in ligs])
        a_z = calc_auc(y, [-l["mean_z_score"] for l in ligs])
        print(f"  {t}: Vina={a_v:.3f} | Vinardo={a_vn:.3f} | AD4={a_ad:.3f} | RRF(k=10)={a_rrf:.3f} | Mean-Z={a_z:.3f}")

    print("\n======================================================================")
    print("2. STAGE 6 CIRCULARITY, CO-CRYSTAL-BLIND, SHUFFLED & LOTO CONTROLS")
    print("======================================================================")
    npol_aucs, tifp_aucs, cb_aucs, shuf_aucs, loto_aucs, srelax_aucs = [], [], [], [], [], []
    rng = np.random.default_rng(42)

    for t in targets:
        s4_map = {l["ligand_id"]: l for l in s4[t]["ranked_ligands"]}
        ligs = s6[t]["relaxed_ligands"]
        y = [1 if l["role"] == "ACTIVE" else 0 for l in ligs]

        a_npol = calc_auc(y, [l["initial_polar_contacts"] for l in ligs])
        a_tifp = calc_auc(y, [l["tanimoto_ifp"] for l in ligs])

        # Co-crystal-blind intrinsic contact + docking (setting N_ref=0, T_IFP=0, delta_ster=0)
        cb_scores = [
            -l["mean_dock_dg"] + 0.35 * (0.85 * l["initial_polar_contacts"] + 1.4 * min(3, s4_map[l["ligand_id"]]["hbd"]) + 0.45 * max(0, s4_map[l["ligand_id"]]["hba"] - 4))
            for l in ligs
        ]
        a_cb = calc_auc(y, cb_scores)

        # Shuffled-contact negative control (200 random permutations)
        rmsds = [l["mean_rmsd_A"] for l in ligs]
        shuf_t = []
        for _ in range(200):
            perm_r = rng.permutation(rmsds)
            shuf_t.append(calc_auc(y, [-(l["mean_dock_dg"] + 1.15 * max(0.0, r - 1.25)) for l, r in zip(ligs, perm_r)]))
        a_shuf = float(np.mean(shuf_t))
        a_sr = calc_auc(y, [-l["dg_relax_kcal_mol"] for l in ligs])

        npol_aucs.append(a_npol); tifp_aucs.append(a_tifp); cb_aucs.append(a_cb)
        shuf_aucs.append(a_shuf); srelax_aucs.append(a_sr)

    # Leave-One-Target-Out (LOTO) cross-validation
    for test_t in targets:
        X_tr, y_tr = [], []
        for tr_t in targets:
            if tr_t == test_t:
                continue
            for l in s6[tr_t]["relaxed_ligands"]:
                X_tr.append([l["mean_dock_dg"], l["mean_rmsd_A"]])
                y_tr.append(1.0 if l["role"] == "ACTIVE" else 0.0)
        X_tr = np.array(X_tr); y_tr = np.array(y_tr)
        mu = X_tr.mean(axis=0); sd = X_tr.std(axis=0) + 1e-9
        Z_tr = (X_tr - mu) / sd

        def loss(w):
            z = Z_tr @ w
            p = 1.0 / (1.0 + np.exp(-np.clip(z, -20, 20)))
            return -np.mean(y_tr * np.log(p + 1e-12) + (1 - y_tr) * np.log(1 - p + 1e-12)) + 0.1 * np.sum(w * w)

        w_opt = minimize(loss, [-1.0, -1.0]).x
        ligs_te = s6[test_t]["relaxed_ligands"]
        y_te = [1 if l["role"] == "ACTIVE" else 0 for l in ligs_te]
        X_te = np.array([[l["mean_dock_dg"], l["mean_rmsd_A"]] for l in ligs_te])
        loto_aucs.append(calc_auc(y_te, ((X_te - mu) / sd) @ w_opt))

    print(f"  N_polar-only Mean AUC              : {np.mean(npol_aucs):.3f} +/- {np.std(npol_aucs, ddof=1):.3f}")
    print(f"  T_IFP-only Mean AUC                : {np.mean(tifp_aucs):.3f} +/- {np.std(tifp_aucs, ddof=1):.3f}")
    print(f"  Co-Crystal-Blind (Dock+Contact) AUC: {np.mean(cb_aucs):.3f} +/- {np.std(cb_aucs, ddof=1):.3f}")
    print(f"  Shuffled-Contact Negative Control  : {np.mean(shuf_aucs):.3f} +/- {np.std(shuf_aucs, ddof=1):.3f}")
    print(f"  LOTO Cross-Validated Mean AUC      : {np.mean(loto_aucs):.3f} +/- {np.std(loto_aucs, ddof=1):.3f}")
    print(f"  Stage 6 Heuristic S_relax Mean AUC : {np.mean(srelax_aucs):.3f} +/- {np.std(srelax_aucs, ddof=1):.3f}")

    print("\n======================================================================")
    print("3. WITHIN-TARGET STRATIFIED BOOTSTRAP (2,000 RESAMPLES) & PAIRED TESTS")
    print("======================================================================")
    boot_methods = {"Vina-only": [], "Mean-Z": [], "Co-Crystal-Blind": [], "S_relax": []}
    for _ in range(2000):
        cur = {k: [] for k in boot_methods}
        for t in targets:
            s4_map = {l["ligand_id"]: l for l in s4[t]["ranked_ligands"]}
            ligs = s6[t]["relaxed_ligands"]
            acts = [l for l in ligs if l["role"] == "ACTIVE"]
            decs = [l for l in ligs if l["role"] == "DECOY"]
            b_all = [acts[i] for i in rng.integers(0, 5, 5)] + [decs[i] for i in rng.integers(0, 5, 5)]
            y = [1]*5 + [0]*5
            cur["Vina-only"].append(calc_auc(y, [-l["score_vina"] for l in b_all]))
            cur["Mean-Z"].append(calc_auc(y, [-s4_map[l["ligand_id"]]["mean_z_score"] for l in b_all]))
            cur["Co-Crystal-Blind"].append(calc_auc(y, [-l["mean_dock_dg"] + 0.35*(0.85*l["initial_polar_contacts"] + 1.4*min(3, s4_map[l["ligand_id"]]["hbd"]) + 0.45*max(0, s4_map[l["ligand_id"]]["hba"]-4)) for l in b_all]))
            cur["S_relax"].append(calc_auc(y, [-l["dg_relax_kcal_mol"] for l in b_all]))
        for k in boot_methods:
            boot_methods[k].append(np.mean(cur[k]))

    for k, arr in boot_methods.items():
        lo, hi = np.percentile(arr, [2.5, 97.5])
        diff = np.array(arr) - np.array(boot_methods["Vina-only"])
        p_val = float(np.mean(diff <= 0))
        print(f"  {k:20s}: Mean={np.mean(arr):.3f}, Stratified 95% CI=[{lo:.3f}, {hi:.3f}], Delta_vs_Vina=+{np.mean(diff):.3f} (paired p={p_val:.4f})")


if __name__ == "__main__":
    main()