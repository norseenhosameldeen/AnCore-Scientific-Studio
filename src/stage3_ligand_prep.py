#!/usr/bin/env python3
"""
AnCore-Scientific-Studio v1.2 — Stage 3: Multi-Molecule 3D Ligand Library Preparation
(RDKit Salt Stripping, Uncharging, ETKDGv3 3D Conformer Embedding, 500-Step MMFF94s
Force-Field Relaxation, Physicochemical & Bickerton QED Profiling, and Torsion-Tree PDBQT Compilation).

Save location in your local project:
    AnCore_v1.2/src/stage3_ligand_prep.py

CLI usage:
  1) Run all 50 benchmark ligands across all 5 targets (1IEP, 3ERT, 1UYG, 1M17, 1STP) (DEFAULT):
     python3 src/stage3_ligand_prep.py --run-all

  2) Include the 7-compound 4WKQ custom target cohort (57 total molecules):
     python3 src/stage3_ligand_prep.py --run-all --include-4wkq

  3) Run a single target (or upload a custom multi-molecule .sdf / .smi library):
     python3 src/stage3_ligand_prep.py --pdb-id 1IEP
     python3 src/stage3_ligand_prep.py --pdb-id 4WKQ --sdf my_compounds.sdf
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
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, QED, rdMolDescriptors

# Cross-version RDKit MolStandardize import (Python 3.10–3.12 / RDKit 2022–2025 compatible)
try:
    from rdkit.Chem.MolStandardize import rdMolStandardize
except ImportError:
    try:
        from rdkit.Chem import rdMolStandardize  # type: ignore
    except ImportError:
        rdMolStandardize = None

try:
    _THIS_DIR = Path(__file__).resolve().parent
except NameError:
    _THIS_DIR = Path.cwd() / "src"
for _p in (_THIS_DIR, Path.cwd(), Path.cwd() / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

BENCHMARK_IDS = ["1IEP", "3ERT", "1UYG", "1M17", "1STP"]

# ── Built-In 50-Compound Benchmark Library (5 Actives + 5 Property-Matched ZINC Decoys per Target)
#    Plus the 7-Compound De Novo 4WKQ Custom Cohort ──────────────────────────────────────────────────
BENCHMARK_LIBRARIES: Dict[str, List[dict]] = {
    "1IEP": [
        {"id": "ABL_ACT_01", "name": "Imatinib (ABL_ACT_01)",      "role": "ACTIVE", "potency": "IC50 = 13 nM",  "smiles": "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc4nccc(n4)c5cccnc5"},
        {"id": "ABL_ACT_02", "name": "Nilotinib (ABL_ACT_02)",     "role": "ACTIVE", "potency": "IC50 = 15 nM",  "smiles": "Cc1cn(cn1)c2cc(NC(=O)c3ccc(C)c(Nc4nccc(n4)c5cccnc5)c3)cc(c2)C(F)(F)F"},
        {"id": "ABL_ACT_03", "name": "Ponatinib (ABL_ACT_03)",     "role": "ACTIVE", "potency": "IC50 = 0.4 nM", "smiles": "Cc1ccc(C(=O)Nc2ccc(CN3CCN(C)CC3)c(c2)C(F)(F)F)cc1C#Cc4cnc5cccnn45"},
        {"id": "ABL_ACT_04", "name": "Dasatinib (ABL_ACT_04)",     "role": "ACTIVE", "potency": "IC50 = 0.8 nM", "smiles": "Cc1cccc(Cl)c1NC(=O)c2cnc(Nc3cc(nc(C)n3)N4CCN(CCO)CC4)s2"},
        {"id": "ABL_ACT_05", "name": "Bosutinib (ABL_ACT_05)",     "role": "ACTIVE", "potency": "IC50 = 1.2 nM", "smiles": "COc1cc2c(Nc3ccc(Cl)c(Cl)c3OC)c(C#N)cnc2cc1OCCCN4CCN(C)CC4"},
        {"id": "ABL_DEC_01", "name": "ZINC_Decoy_A1",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "COc1ccc(cc1)S(=O)(=O)N2CCN(CC2)C(=O)c3ccc4ccccc4c3"},
        {"id": "ABL_DEC_02", "name": "ZINC_Decoy_A2",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "Cc1ccc(cc1)C(=O)N2CCC(CC2)CNC(=O)c3cc4ccccc4o3"},
        {"id": "ABL_DEC_03", "name": "ZINC_Decoy_A3",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "COc1ccc(cc1OC)C(=O)N2CCN(CC2)c3nc4ccccc4s3"},
        {"id": "ABL_DEC_04", "name": "ZINC_Decoy_A4",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "Cc1cc(C)nc(SCC(=O)Nc2ccc(cc2)S(=O)(=O)N3CCOCC3)n1"},
        {"id": "ABL_DEC_05", "name": "ZINC_Decoy_A5",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "Fc1ccc(cc1)C(=O)N2CCN(CC2)C(=O)Cc3ccc4OCOc4c3"},
    ],
    "3ERT": [
        {"id": "ERA_ACT_01", "name": "4-Hydroxytamoxifen (OHT)",   "role": "ACTIVE", "potency": "Ki = 0.5 nM",   "smiles": "CC/C(=C(\\c1ccccc1)/c2ccc(OCCN(C)C)cc2)/c3ccc(O)cc3"},
        {"id": "ERA_ACT_02", "name": "Raloxifene (ERA_ACT_02)",    "role": "ACTIVE", "potency": "Ki = 0.4 nM",   "smiles": "Oc1ccc2c(c1)sc(c3ccc(O)cc3)c2C(=O)c4ccc(OCCN5CCCCC5)cc4"},
        {"id": "ERA_ACT_03", "name": "Lasofoxifene (ERA_ACT_03)",  "role": "ACTIVE", "potency": "Ki = 1.5 nM",   "smiles": "Oc1ccc2CCC(c3ccccc3)C(c4ccc(OCCN5CCCC5)cc4)c2c1"},
        {"id": "ERA_ACT_04", "name": "Bazedoxifene (ERA_ACT_04)",  "role": "ACTIVE", "potency": "Ki = 0.1 nM",   "smiles": "Cc1c(c2ccc(O)cc2)n(Cc3ccc(OCCN4CCCCCC4)cc3)c5ccc(O)cc15"},
        {"id": "ERA_ACT_05", "name": "17beta-Estradiol",           "role": "ACTIVE", "potency": "Kd = 0.1 nM",   "smiles": "CC12CCC3c4ccc(O)cc4CCC3C1CCC2O"},
        {"id": "ERA_DEC_01", "name": "ZINC_Decoy_E1",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "COc1ccc(cc1)C2CCN(CC2)C(=O)c3ccc(cc3)S(=O)(=O)N(C)C"},
        {"id": "ERA_DEC_02", "name": "ZINC_Decoy_E2",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "Cc1ccc(cc1)S(=O)(=O)N2CCCC2C(=O)Nc3ccc(Cl)cc3"},
        {"id": "ERA_DEC_03", "name": "ZINC_Decoy_E3",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "COc1ccc(cc1)NC(=O)C2CCCN(C2)S(=O)(=O)c3ccc(F)cc3"},
        {"id": "ERA_DEC_04", "name": "ZINC_Decoy_E4",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "CCOC(=O)c1ccc(NC(=O)CN2CCN(CC2)c3ccccc3)cc1"},
        {"id": "ERA_DEC_05", "name": "ZINC_Decoy_E5",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "Cc1ccc(NC(=O)c2CCCN(C2)C(=O)c3ccco3)cc1C"},
    ],
    "1UYG": [
        {"id": "HSP_ACT_01", "name": "PU3 purine (HSP_ACT_01)",    "role": "ACTIVE", "potency": "IC50 = 15 nM",  "smiles": "COc1cc(OC)c(Cl)c(c1)Cc2nc(N)nc3n(CCCC#C)cnc23"},
        {"id": "HSP_ACT_02", "name": "BIIB021 (HSP_ACT_02)",       "role": "ACTIVE", "potency": "IC50 = 5.1 nM", "smiles": "COc1nc(C)c(Cc2cnc3nc(N)ncn23)c(Cl)c1C"},
        {"id": "HSP_ACT_03", "name": "Ganetespib (HSP_ACT_03)",    "role": "ACTIVE", "potency": "IC50 = 4.0 nM", "smiles": "CC(C)c1cc(c(O)cc1O)c2nn(c3ccc4ccccc4c3)c(=O)[nH]2"},
        {"id": "HSP_ACT_04", "name": "Luminespib (HSP_ACT_04)",    "role": "ACTIVE", "potency": "IC50 = 1.7 nM", "smiles": "CCNC(=O)c1noc(c2cc(C(C)C)c(O)cc2O)c1c3ccc(CN4CCOCC4)cc3"},
        {"id": "HSP_ACT_05", "name": "Onalespib (HSP_ACT_05)",     "role": "ACTIVE", "potency": "IC50 = 0.7 nM", "smiles": "CC(C)c1cc(C(=O)N2Cc3ccc(CN4CCN(C)CC4)cc3C2)c(O)cc1O"},
        {"id": "HSP_DEC_01", "name": "ZINC_Decoy_H1",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "COc1ccc(cc1)C(=O)Nc2cccc(c2)S(=O)(=O)N3CCCC3"},
        {"id": "HSP_DEC_02", "name": "ZINC_Decoy_H2",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "Cc1ccc(cc1)S(=O)(=O)NCC(=O)N2CCN(CC2)c3ccccc3"},
        {"id": "HSP_DEC_03", "name": "ZINC_Decoy_H3",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "CCOC(=O)c1ccc(cc1)N2CCN(CC2)C(=O)c3ccc(F)cc3"},
        {"id": "HSP_DEC_04", "name": "ZINC_Decoy_H4",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "COc1ccc(cc1)NC(=O)c2ccc(cc2)S(=O)(=O)N3CCOCC3"},
        {"id": "HSP_DEC_05", "name": "ZINC_Decoy_H5",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "Cc1cccc(c1)N2CCN(CC2)C(=O)CNS(=O)(=O)c3ccc(Cl)cc3"},
    ],
    "1M17": [
        {"id": "EGF_ACT_01", "name": "Erlotinib (EGF_ACT_01)",     "role": "ACTIVE", "potency": "IC50 = 2.0 nM",  "smiles": "COCCOc1cc2ncnc(Nc3cccc(c3)C#C)c2cc1OCCOC"},
        {"id": "EGF_ACT_02", "name": "Gefitinib (EGF_ACT_02)",     "role": "ACTIVE", "potency": "IC50 = 33 nM",   "smiles": "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN4CCOCC4"},
        {"id": "EGF_ACT_03", "name": "Afatinib (EGF_ACT_03)",      "role": "ACTIVE", "potency": "IC50 = 0.5 nM",  "smiles": "CN(C)C/C=C/C(=O)Nc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cc1OC4CCOC4"},
        {"id": "EGF_ACT_04", "name": "Lapatinib (EGF_ACT_04)",     "role": "ACTIVE", "potency": "IC50 = 10.8 nM", "smiles": "CS(=O)(=O)CCNCc1ccc(o1)c2ccc3ncnc(Nc4ccc(OCc5cccc(F)c5)c(Cl)c4)c3c2"},
        {"id": "EGF_ACT_05", "name": "Osimertinib (EGF_ACT_05)",   "role": "ACTIVE", "potency": "IC50 = 12 nM",   "smiles": "COc1cc(N(C)CCN(C)C)c(NC(=O)C=C)cc1Nc2nccc(n2)c3cn(C)c4ccccc34"},
        {"id": "EGF_DEC_01", "name": "ZINC_Decoy_G1",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "COc1ccc(cc1)C(=O)N2CCN(CC2)S(=O)(=O)c3ccc(C)cc3"},
        {"id": "EGF_DEC_02", "name": "ZINC_Decoy_G2",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "Cc1ccc(cc1)NC(=O)CN2CCN(CC2)C(=O)c3ccc(Cl)cc3"},
        {"id": "EGF_DEC_03", "name": "ZINC_Decoy_G3",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "CCOC(=O)c1ccc(NC(=O)c2ccc(cc2)S(=O)(=O)N3CCCC3)cc1"},
        {"id": "EGF_DEC_04", "name": "ZINC_Decoy_G4",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "COc1ccc(cc1OC)S(=O)(=O)N2CCN(CC2)C(=O)c3ccco3"},
        {"id": "EGF_DEC_05", "name": "ZINC_Decoy_G5",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "Fc1ccc(cc1)S(=O)(=O)N2CCN(CC2)C(=O)Nc3ccc(C)cc3"},
    ],
    "1STP": [
        {"id": "STP_ACT_01", "name": "D-Biotin (STP_ACT_01)",      "role": "ACTIVE", "potency": "Kd = 0.04 pM",   "smiles": "O=C(O)CCCCC1SCC2NC(=O)NC12"},
        {"id": "STP_ACT_02", "name": "Desthiobiotin (STP_ACT_02)", "role": "ACTIVE", "potency": "Kd = 1.0 nM",    "smiles": "CC1NC(=O)NC1CCCCCC(=O)O"},
        {"id": "STP_ACT_03", "name": "2-Iminobiotin (STP_ACT_03)", "role": "ACTIVE", "potency": "Kd = 3.5 nM",    "smiles": "NC1=NC2CSCC(CCCCC(=O)O)C2N1"},
        {"id": "STP_ACT_04", "name": "Diaminobiotin (STP_ACT_04)", "role": "ACTIVE", "potency": "Kd = 12 nM",     "smiles": "NC1C(CCCCC(=O)O)SC2CN12"},
        {"id": "STP_ACT_05", "name": "alpha-Lipoic acid (STP_ACT_05)", "role": "ACTIVE", "potency": "Kd = 700 nM", "smiles": "O=C(O)CCCCC1CCSS1"},
        {"id": "STP_DEC_01", "name": "ZINC_Decoy_S1",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "Cc1ccc(cc1)S(=O)(=O)NC2CCCC2=O"},
        {"id": "STP_DEC_02", "name": "ZINC_Decoy_S2",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "O=C1CCC(C(=O)N1c2ccc(Cl)cc2)c3ccccc3"},
        {"id": "STP_DEC_03", "name": "ZINC_Decoy_S3",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "O=C1CC(C(=O)N1c2ccc(F)cc2)c3ccc(F)cc3"},
        {"id": "STP_DEC_04", "name": "ZINC_Decoy_S4",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "COc1ccc(cc1)N2C(=O)CCC2C(=O)O"},
        {"id": "STP_DEC_05", "name": "ZINC_Decoy_S5",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "Cc1ccc(NC(=O)C2CCCN2C(=O)C)cc1"},
    ],
    "4WKQ": [
        {"id": "WKQ_ACT_01", "name": "Gefitinib (WKQ_ACT_01)",     "role": "ACTIVE", "potency": "IC50 = 33 nM",   "smiles": "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN4CCOCC4"},
        {"id": "WKQ_ACT_02", "name": "Erlotinib (WKQ_ACT_02)",     "role": "ACTIVE", "potency": "IC50 = 2.0 nM",  "smiles": "COCCOc1cc2ncnc(Nc3cccc(c3)C#C)c2cc1OCCOC"},
        {"id": "WKQ_ACT_03", "name": "Osimertinib (WKQ_ACT_03)",   "role": "ACTIVE", "potency": "IC50 = 12 nM",   "smiles": "COc1cc(N(C)CCN(C)C)c(NC(=O)C=C)cc1Nc2nccc(n2)c3cn(C)c4ccccc34"},
        {"id": "WKQ_ACT_04", "name": "Imatinib (WKQ_ACT_04)",      "role": "ACTIVE", "potency": "IC50 = 13 nM",   "smiles": "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc4nccc(n4)c5cccnc5"},
        {"id": "WKQ_DEC_01", "name": "ZINC_Decoy_W1",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "COc1ccc(cc1)C(=O)N2CCN(CC2)S(=O)(=O)c3ccc(C)cc3"},
        {"id": "WKQ_DEC_02", "name": "ZINC_Decoy_W2",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "Cc1ccc(cc1)NC(=O)CN2CCN(CC2)C(=O)c3ccc(Cl)cc3"},
        {"id": "WKQ_DEC_03", "name": "ZINC_Decoy_W3",              "role": "DECOY",  "potency": "Decoy (>10 uM)", "smiles": "Fc1ccc(cc1)S(=O)(=O)N2CCN(CC2)C(=O)Nc3ccc(C)cc3"},
    ],
}


def standardize_and_strip_salts(raw_mol: Chem.Mol) -> Chem.Mol:
    """
    Applies RDKit LargestFragmentChooser + Uncharger when `rdMolStandardize` is available,
    falling back to `Chem.GetMolFrags(asMols=True)` largest heavy-atom fragment selection.
    """
    if rdMolStandardize is not None:
        clean_mol = rdMolStandardize.LargestFragmentChooser().choose(raw_mol)
        clean_mol = rdMolStandardize.Uncharger().uncharge(clean_mol)
        return clean_mol
    frags = Chem.GetMolFrags(raw_mol, asMols=True, sanitizeFrags=True)
    return max(frags, key=lambda m: m.GetNumHeavyAtoms()) if frags else raw_mol


def strip_pdbqt_model_tags(pdbqt_text: str) -> str:
    """Strips MODEL / ENDMDL tags so AutoDock Vina's C++ parser accepts single-ligand PDBQT files."""
    return "".join(
        line + "\n" for line in pdbqt_text.splitlines() if not line.startswith(("MODEL", "ENDMDL"))
    )


def translate_mol_to_center(mol: Chem.Mol, target_center: List[float]) -> None:
    """Translates an RDKit conformer centroid to the receptor binding-pocket center (cx, cy, cz)."""
    conf = mol.GetConformer()
    pts = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())], dtype=float)
    centroid = pts.mean(axis=0)
    shift = np.array(target_center, dtype=float) - centroid
    for i in range(mol.GetNumAtoms()):
        p = pts[i] + shift
        conf.SetAtomPosition(i, (float(p[0]), float(p[1]), float(p[2])))


def compile_ligand_to_pdbqt(mol_h: Chem.Mol, temp_sdf: Path, out_pdbqt: Path) -> str:
    """Compiles a 3D protonated RDKit molecule into a flexible torsion-tree .pdbqt file."""
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy

        prep = MoleculePreparation()
        setup_list = prep.prepare(mol_h)
        pdbqt_str, is_ok, _ = PDBQTWriterLegacy.write_string(setup_list[0])
        if is_ok and ("ATOM" in pdbqt_str or "HETATM" in pdbqt_str):
            out_pdbqt.write_text(strip_pdbqt_model_tags(pdbqt_str), encoding="utf-8")
            return "Meeko Torsion-Tree Compiler"
    except Exception:
        pass

    obabel_bin = shutil.which("obabel")
    if obabel_bin:
        cmd = [obabel_bin, "-isdf", str(temp_sdf), "-opdbqt", "-O", str(out_pdbqt)]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode == 0 and out_pdbqt.is_file() and out_pdbqt.stat().st_size > 100:
            clean_txt = strip_pdbqt_model_tags(out_pdbqt.read_text(encoding="utf-8", errors="ignore"))
            out_pdbqt.write_text(clean_txt, encoding="utf-8")
            return "OpenBabel Torsion-Tree Compiler"

    raise RuntimeError(f"Failed to compile {out_pdbqt.name}: ensure `meeko` or `obabel` is installed.")


def load_custom_sdf_or_smiles(sdf_path: Optional[Path], smiles_file: Optional[Path], pdb_id: str) -> List[dict]:
    """Loads candidate molecules from a multi-molecule .sdf file, .smi file, or built-in benchmark cohort."""
    if sdf_path is not None and sdf_path.is_file():
        suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
        records = []
        for idx, mol in enumerate(suppl, start=1):
            if mol is None:
                continue
            name = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else ""
            for tag in ("NAME", "TITLE", "ID", "COMPOUND_NAME", "ZINC_ID", "CHEMBL_ID"):
                if not name and mol.HasProp(tag):
                    name = mol.GetProp(tag).strip()
            if not name:
                name = f"{pdb_id}_LIG_{idx:02d}"
            role = mol.GetProp("ROLE").strip().upper() if mol.HasProp("ROLE") else ("DECOY" if "DECOY" in name.upper() else "ACTIVE")
            potency = mol.GetProp("POTENCY").strip() if mol.HasProp("POTENCY") else "Custom Upload"
            records.append({
                "id": f"{pdb_id}_{idx:02d}",
                "name": name,
                "role": role,
                "potency": potency,
                "smiles": Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=True),
            })
        if records:
            return records

    if smiles_file is not None and smiles_file.is_file():
        records = []
        for idx, ln in enumerate(smiles_file.read_text(encoding="utf-8").splitlines(), start=1):
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            parts = ln.split()
            smi = parts[0]
            name = parts[1] if len(parts) > 1 else f"{pdb_id}_LIG_{idx:02d}"
            role = "DECOY" if "DECOY" in name.upper() else "ACTIVE"
            records.append({
                "id": f"{pdb_id}_{idx:02d}",
                "name": name,
                "role": role,
                "potency": "Custom SMILES",
                "smiles": smi,
            })
        if records:
            return records

    return list(BENCHMARK_LIBRARIES.get(pdb_id.upper(), BENCHMARK_LIBRARIES["1IEP"]))


def prepare_single_ligand(
    entry: dict,
    pocket_center: List[float],
    pdbqt_dir: Path,
    seed: int = 42,
) -> Tuple[dict, Chem.Mol]:
    """
    Standardizes SMILES, strips salts, neutralizes/enumerates protomer state, embeds 3D coordinates
    via ETKDGv3, relaxes for 500 steps under MMFF94s (Eq. 5), computes Bickerton QED (Eq. 6) and
    Lipinski descriptors, centers in the binding pocket, and compiles a torsion-tree .pdbqt.
    """
    lig_id = entry["id"]
    raw_mol = Chem.MolFromSmiles(entry["smiles"])
    if raw_mol is None:
        raise ValueError(f"Invalid SMILES for {lig_id}: {entry['smiles']}")

    # 1. Salt stripping & canonical neutralization
    clean_mol = standardize_and_strip_salts(raw_mol)
    can_smiles = Chem.MolToSmiles(clean_mol, canonical=True)

    # 2. Physicochemical descriptors & Bickerton QED (Eq. 6) on standardized parent molecule
    mw = float(Descriptors.MolWt(clean_mol))
    clogp = float(Descriptors.MolLogP(clean_mol))
    tpsa = float(Descriptors.TPSA(clean_mol))
    hbd = int(Descriptors.NumHDonors(clean_mol))
    hba = int(Descriptors.NumHAcceptors(clean_mol))
    rotb = int(Descriptors.NumRotatableBonds(clean_mol))
    n_heavy = int(clean_mol.GetNumHeavyAtoms())
    qed_score = float(QED.qed(clean_mol))
    formula = rdMolDescriptors.CalcMolFormula(clean_mol)
    ro5_violations = int(mw > 500.0) + int(clogp > 5.0) + int(hbd > 5) + int(hba > 10)

    # 3. 3D ETKDGv3 distance-geometry embedding + 500-step MMFF94s force-field relaxation (Eq. 5)
    mol_h = Chem.AddHs(clean_mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    res_emb = AllChem.EmbedMolecule(mol_h, params)
    if res_emb != 0:
        params2 = AllChem.ETKDGv2()
        params2.randomSeed = int(seed)
        if AllChem.EmbedMolecule(mol_h, params2) != 0:
            AllChem.EmbedMolecule(mol_h, useRandomCoords=True, randomSeed=int(seed))

    forcefield_used = "MMFF94s"
    mmff_energy = 0.0
    if AllChem.MMFFHasAllMoleculeParams(mol_h):
        AllChem.MMFFOptimizeMolecule(mol_h, mmffVariant="MMFF94s", maxIters=500)
        mp = AllChem.MMFFGetMoleculeProperties(mol_h, mmffVariant="MMFF94s")
        if mp is not None:
            ff = AllChem.MMFFGetMoleculeForceField(mol_h, mp)
            if ff is not None:
                mmff_energy = float(ff.CalcEnergy())
    else:
        forcefield_used = "UFF"
        AllChem.UFFOptimizeMolecule(mol_h, maxIters=500)
        ff = AllChem.UFFGetMoleculeForceField(mol_h)
        if ff is not None:
            mmff_energy = float(ff.CalcEnergy())

    # 4. Translate conformer centroid to the Stage 1 pocket center (cx, cy, cz)
    translate_mol_to_center(mol_h, pocket_center)

    mol_h.SetProp("_Name", entry["name"])
    mol_h.SetProp("LIGAND_ID", lig_id)
    mol_h.SetProp("ROLE", entry["role"])
    mol_h.SetProp("CANONICAL_SMILES", can_smiles)
    mol_h.SetProp("FORMULA", formula)
    mol_h.SetProp("MW", f"{mw:.2f}")
    mol_h.SetProp("CLOGP", f"{clogp:.2f}")
    mol_h.SetProp("TPSA", f"{tpsa:.2f}")
    mol_h.SetProp("HBD", str(hbd))
    mol_h.SetProp("HBA", str(hba))
    mol_h.SetProp("ROTB", str(rotb))
    mol_h.SetProp("QED", f"{qed_score:.3f}")
    mol_h.SetProp("MMFF94S_ENERGY_KCAL", f"{mmff_energy:.2f}")
    mol_h.SetProp("EXPERIMENTAL_POTENCY", entry.get("potency", "N/A"))

    # 5. Export individual .sdf and compile torsion-tree .pdbqt
    temp_sdf = pdbqt_dir / f"{lig_id}_prepared.sdf"
    with Chem.SDWriter(str(temp_sdf)) as w:
        w.write(mol_h)

    out_pdbqt = pdbqt_dir / f"{lig_id}.pdbqt"
    compiler_name = compile_ligand_to_pdbqt(mol_h, temp_sdf, out_pdbqt)

    rec_dict = {
        "ligand_id": lig_id,
        "name": entry["name"],
        "role": entry["role"],
        "canonical_smiles": can_smiles,
        "formula": formula,
        "heavy_atoms": n_heavy,
        "mw": round(mw, 1),
        "clogp": round(clogp, 2),
        "tpsa_A2": round(tpsa, 1),
        "hbd": hbd,
        "hba": hba,
        "hbd_hba": f"{hbd}/{hba}",
        "rot_bonds": rotb,
        "qed": round(qed_score, 3),
        "lipinski_violations": ro5_violations,
        "lipinski_pass": ro5_violations <= 1,
        "forcefield": forcefield_used,
        "mmff94s_energy_kcal_mol": round(mmff_energy, 2),
        "experimental_potency": entry.get("potency", "N/A"),
        "pdbqt_compiler": compiler_name,
        "pdbqt_path": str(out_pdbqt),
        "sdf_path": str(temp_sdf),
    }
    return rec_dict, mol_h


def run_stage3_for_target(
    pdb_id: str,
    workspace_dir: Path,
    sdf_path: Optional[Path] = None,
    smiles_file: Optional[Path] = None,
    seed: int = 42,
) -> dict:
    """Runs Stage 3 3D conformer embedding, MMFF94s relaxation, QED profiling, and PDBQT compilation."""
    workspace_dir = workspace_dir.resolve()
    s1_manifest_path = workspace_dir / "stage1_manifest.json"
    if not s1_manifest_path.is_file():
        from stage1_target_prep import run_stage1_target_prep
        run_stage1_target_prep(pdb_id=pdb_id, out_dir=workspace_dir)

    s1 = json.loads(s1_manifest_path.read_text(encoding="utf-8"))
    pocket_center = [
        s1["docking_box"]["center_x"],
        s1["docking_box"]["center_y"],
        s1["docking_box"]["center_z"],
    ]

    s3_dir = workspace_dir / "03_Stage3_Ligand_Library"
    pdbqt_dir = s3_dir / "meeko_pdbqt_inputs"
    pdbqt_dir.mkdir(parents=True, exist_ok=True)

    raw_entries = load_custom_sdf_or_smiles(sdf_path, smiles_file, pdb_id.upper())
    prepared_records = []
    mol_objects = []

    for idx, entry in enumerate(raw_entries):
        rec, mol_h = prepare_single_ligand(entry, pocket_center, pdbqt_dir, seed=seed + idx)
        prepared_records.append(rec)
        mol_objects.append(mol_h)

    combined_sdf = s3_dir / "standardized_ligands_3D_MMFF94s.sdf"
    with Chem.SDWriter(str(combined_sdf)) as writer:
        for m in mol_objects:
            writer.write(m)

    csv_path = s3_dir / "ligand_physicochemical_properties.csv"
    csv_fields = [
        "ligand_id", "name", "role", "mw", "clogp", "tpsa_A2", "hbd_hba",
        "rot_bonds", "qed", "mmff94s_energy_kcal_mol", "lipinski_pass", "experimental_potency",
    ]
    with open(csv_path, "w", newline="", encoding="utf-8") as fc:
        w = csv.DictWriter(fc, fieldnames=csv_fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(prepared_records)

    actives = [r for r in prepared_records if r["role"] == "ACTIVE"]
    decoys = [r for r in prepared_records if r["role"] == "DECOY"]

    manifest = {
        "stage": 3,
        "status": "COMPLETED",
        "pdb_id": pdb_id.upper(),
        "total_ligands": len(prepared_records),
        "active_count": len(actives),
        "decoy_count": len(decoys),
        "mean_mw": round(float(np.mean([r["mw"] for r in prepared_records])), 1),
        "mean_clogp": round(float(np.mean([r["clogp"] for r in prepared_records])), 2),
        "mean_tpsa_A2": round(float(np.mean([r["tpsa_A2"] for r in prepared_records])), 1),
        "mean_rot_bonds": round(float(np.mean([r["rot_bonds"] for r in prepared_records])), 1),
        "active_mean_rot_bonds": round(float(np.mean([r["rot_bonds"] for r in actives])), 1) if actives else 0.0,
        "decoy_mean_rot_bonds": round(float(np.mean([r["rot_bonds"] for r in decoys])), 1) if decoys else 0.0,
        "mean_qed": round(float(np.mean([r["qed"] for r in prepared_records])), 3),
        "mean_mmff94s_energy_kcal_mol": round(float(np.mean([r["mmff94s_energy_kcal_mol"] for r in prepared_records])), 2),
        "ligands": prepared_records,
        "files": {
            "standardized_sdf": str(combined_sdf),
            "physicochemical_csv": str(csv_path),
            "pdbqt_dir": str(pdbqt_dir),
        },
    }

    (workspace_dir / "stage3_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="AnCore v1.2 Stage 3: Multi-Molecule 3D Ligand Library Preparation")
    parser.add_argument("--pdb-id", type=str, default=None, help="Single 4-character PDB ID (omit to run all 5 benchmark targets)")
    parser.add_argument("--run-all", action="store_true", help="Run all 50 ligands across the 5 benchmark targets")
    parser.add_argument("--include-4wkq", action="store_true", help="Also prepare the 7-compound 4WKQ library")
    parser.add_argument("--sdf", type=str, default=None, help="Optional path to a custom multi-molecule .sdf file")
    parser.add_argument("--smiles", type=str, default=None, help="Optional path to a custom .smi SMILES file")
    parser.add_argument("--out-dir", type=str, default="workspaces", help="Base workspaces directory (default: workspaces)")
    args, _ = parser.parse_known_args()

    base_out = Path(args.out_dir)

    if args.run_all or args.pdb_id is None:
        targets = list(BENCHMARK_IDS)
        if args.include_4wkq:
            targets.append("4WKQ")

        print("=" * 112)
        print("ANCORE v1.2 — STAGE 3 MULTI-TARGET 3D LIGAND CONFORMER PREPARATION (ETKDGv3 + MMFF94s + QED)")
        print("=" * 112)
        print(f"{'PDB':<6} {'Cohort':<12} {'Mean MW':<10} {'cLogP':<8} {'TPSA (Å²)':<11} {'Act/Dec RotB':<15} {'Mean QED':<10} {'Mean E_MMFF94s':<16} {'Status'}")
        print("-" * 112)

        suite_results = []
        all_ligs = []
        for pid in targets:
            m = run_stage3_for_target(pid, base_out / pid)
            suite_results.append(m)
            all_ligs.extend(m["ligands"])
            cohort_str = f"{m['active_count']}A + {m['decoy_count']}D"
            rot_str = f"{m['active_mean_rot_bonds']:.1f} / {m['decoy_mean_rot_bonds']:.1f}"
            print(
                f"{pid:<6} {cohort_str:<12} {m['mean_mw']:<10.1f} {m['mean_clogp']:<8.2f} "
                f"{m['mean_tpsa_A2']:<11.1f} {rot_str:<15} {m['mean_qed']:<10.3f} "
                f"{m['mean_mmff94s_energy_kcal_mol']:<16.2f} 100% READY ({m['total_ligands']}/{m['total_ligands']})"
            )

        print("-" * 112)
        tot_n = len(all_ligs)
        tot_act = sum(1 for r in all_ligs if r["role"] == "ACTIVE")
        tot_dec = sum(1 for r in all_ligs if r["role"] == "DECOY")
        print(
            f"{'SUITE':<6} {f'{tot_act}A + {tot_dec}D':<12} "
            f"{float(np.mean([r['mw'] for r in all_ligs])):<10.1f} "
            f"{float(np.mean([r['clogp'] for r in all_ligs])):<8.2f} "
            f"{float(np.mean([r['tpsa_A2'] for r in all_ligs])):<11.1f} "
            f"{float(np.mean([r['rot_bonds'] for r in all_ligs])):<15.1f} "
            f"{float(np.mean([r['qed'] for r in all_ligs])):<10.3f} "
            f"{float(np.mean([r['mmff94s_energy_kcal_mol'] for r in all_ligs])):<16.2f} "
            f"100% READY ({tot_n}/{tot_n})"
        )
        print("=" * 112)

        suite_json = base_out / "stage3_suite_summary.json"
        suite_json.write_text(json.dumps(suite_results, indent=2), encoding="utf-8")
        print(f"Saved multi-target Stage 3 summary -> {suite_json}")
        return

    sdf_p = Path(args.sdf) if args.sdf else None
    smi_p = Path(args.smiles) if args.smiles else None
    m = run_stage3_for_target(args.pdb_id.upper(), base_out / args.pdb_id.upper(), sdf_path=sdf_p, smiles_file=smi_p)
    print(json.dumps(m, indent=2))


if __name__ == "__main__":
    main()