# AnCore Scientific Studio (v1.2.0)

**A Six-Stage Structure-Based Virtual Screening Workbench Integrating Two-Phase Redocking, Reciprocal Rank Fusion Consensus Scoring, and Three-Zone 300 K Langevin Dynamics Relaxation**

**Authors:** Norseen H. Abdelhakim and Ahmed M. Senary  
**Affiliation:** Department of Medicinal Chemistry, Faculty of Pharmacy, Galala University, New Galala 43713, Egypt  
**Repository:** [github.com/norseenhosameldeen/AnCore-Scientific-Studio](https://github.com/norseenhosameldeen/AnCore-Scientific-Studio)  
**Preprint DOI:** [10.26434/chemrxiv.15012010/v1](https://www.google.com/search?q=https://doi.org/10.26434/chemrxiv.15012010/v1)  
**License:** MIT License

## Overview

**AnCore Scientific Studio (v1.2.0)** is a native macOS desktop workbench and command-line chemoinformatics pipeline for reproducible, multi-ligand structure-based virtual screening (SBVS) and physics-based post-docking pose validation.

Desktop virtual screening tools frequently encounter two practical bottlenecks:
1. **Graphical IPC Serialization Bottleneck:** Python-wrapped webview shells (`pywebview` / Electron bridges) serialize molecular coordinates as JSON strings over inter-process pipes, failing above a ~2 MB payload ceiling when rendering full protein receptors alongside multi-pose ligand ensembles or molecular dynamics (MD) trajectories.
2. **Static Docking False Positives:** Stopping a virtual screen at a single static docking score ignores thermal motion at 300 K, desolvation penalties, and conformational entropy, allowing rigid decoys with opportunistic surface contacts to outrank genuine bioactive inhibitors.

AnCore v1.2.0 resolves both bottlenecks by decoupling a compiled [Rust](https://www.rust-lang.org/) and [Tauri](https://tauri.app/) desktop backend (`ancore-tauri/src-tauri/`) with a hardware-accelerated native macOS `WKWebView` / [3Dmol.js](https://3dmol.csb.pitt.edu/) WebGL viewport (`ancore-tauri/ui/index.html`) from a six-stage Python biophysical engine (`src/`) powered by [PDBFixer](https://github.com/openmm/pdbfixer), [OpenMM](https://openmm.org/), [RDKit](https://www.rdkit.org/), [OpenBabel](https://openbabel.org/), [AutoDock Vina](https://vina.scripps.edu/), and a turnkey Mode C 100-ns explicit-solvent [GROMACS](https://www.gromacs.org/) / OpenMM package exporter.

## Key Features in v1.2.0

* **Zero-Limit Native Rust/Tauri & WebGL Architecture:** Streams receptor structures, multi-pose docking ensembles, and MD trajectory frames directly over an asynchronous localhost HTTP loopback server with zero string-serialization freezes.
* **Stage 1 — Receptor Curation, AMBER14/GBn2 Minimization & Steric Clash Gate:** Automates hetero-atom stripping, [PDBFixer](https://github.com/openmm/pdbfixer) heavy-atom reconstruction, physiological pH 7.4 protonation, 200-step L-BFGS relaxation in [OpenMM](https://openmm.org/) (`amber14-all.xml` + `implicit/gbn2.xml`), non-bonded steric clash gating ($\min d_{ij} \ge 2.10\text{ \AA}$), and 5.0 Å co-crystal pocket bounding-box mapping.
* **Stage 2 — Two-Phase Co-Crystal Redocking Gate ($\text{RMSD} \le 2.00\text{ \AA}$):** Combines Phase A local BFGS energy refinement (`--local_only`) with Phase B focused-box ($\delta_{\text{pad}} = 3.0\text{ \AA}$, `exhaustiveness = 16`) global Monte Carlo search, rescuing high-rotatable-bond kinase ligands from unfocused search failures before permitting library screening.
* **Stage 3 — 3D Ligand Library Embedding & Physicochemical Profiling:** Performs pH 7.4 protomer normalization, [RDKit](https://www.rdkit.org/) ETKDGv3 3D distance-geometry embedding, 500-iteration MMFF94s force-field relaxation, Gasteiger-Marsili partial charge assignment, and Lipinski/QED profiling.
* **Stage 4 — Tri-Scoring Consensus Virtual Screening (RRF + Mean Z-Score):** Evaluates every ligand across three distinct scoring potentials—[AutoDock Vina](https://vina.scripps.edu/) (`Vina`), `Vinardo` (local BFGS rescoring), and `AutoDock4` (`AD4`)—and fuses rankings via Reciprocal Rank Fusion ($\text{RRF}_{k=10}$) and standardized mean Z-score ($\bar{Z}$).
* **Stage 5 — 3D Protein-Ligand Interaction Fingerprinting (IFP) & 2D Pharmacophore Maps:** Computes 3D geometric contacts across hydrogen bonds ($\le 3.5\text{ \AA}$), salt bridges ($\le 4.2\text{ \AA}$), $\pi$-stacking ($\le 4.3\text{ \AA}$), halogen bonds ($\le 3.9\text{ \AA}$), and hydrophobic contacts ($\le 4.0\text{ \AA}$), alongside reference Tanimoto similarity ($T_{\text{IFP}}$) and 300-DPI 2D interaction diagrams.
* **Stage 6 — Three-Zone 300 K OpenMM Langevin Relaxation & Mode C 100-ns Exporter:** Subjects top consensus complexes to intact-topology three-zone 300 K Langevin dynamics ($\gamma = 18.0\text{ ps}^{-1}$, $\Delta t = 2.0\text{ fs}$, 5.0 ns benchmark window) to eject thermally unstable false-positive decoys, and exports ready-to-run 100-ns ($50,000,000$-step) explicit-solvent OpenMM XML and [GROMACS](https://www.gromacs.org/) 2024 simulation bundles (`Mode C`).

## Six-Stage Pipeline Architecture

| Stage | Module Script | Core Scientific Engines | Physical Quality Gate / Output Metrics |
| :--- | :--- | :--- | :--- |
| **Stage 1** | `src/stage1_target_prep.py` | [PDBFixer](https://github.com/openmm/pdbfixer), [OpenMM](https://openmm.org/) (`amber14-all` + `gbn2`), [OpenBabel](https://openbabel.org/) | pH 7.4 protonation; 200-step L-BFGS minimization; steric clash gate $\min d_{ij} \ge 2.10\text{ \AA}$; 5.0 Å pocket box |
| **Stage 2** | `src/stage2_redock_gate.py` | [AutoDock Vina](https://vina.scripps.edu/) 1.2, SciPy Hungarian Kabsch RMSD | Phase A local BFGS + Phase B focused MC ($\delta_{\text{pad}} = 3.0\text{ \AA}$, `exhaustiveness = 16`); gate $\text{RMSD} \le 2.00\text{ \AA}$ |
| **Stage 3** | `src/stage3_ligand_prep.py` | [RDKit](https://www.rdkit.org/) (`ETKDGv3`, `MMFF94s`), [OpenBabel](https://openbabel.org/) | `.SDF`/`.MOL`/SMILES import; 500-iter MMFF94s relaxation; MW, cLogP, TPSA, HBD, HBA, $N_{\text{rot}}$, QED |
| **Stage 4** | `src/stage4_screening.py` | [AutoDock Vina](https://vina.scripps.edu/) (`Vina`, `Vinardo`, `AD4`) | Discovery box ($\delta_{\text{pad}} = 4.0\text{ \AA}$); $\text{RRF}_{k=10}$ consensus score; mean Z-score $\bar{Z}$; $\text{ROC-AUC}$ & $\text{EF}_{20\%}$ |
| **Stage 5** | `src/stage5_interactions.py` | AnCore 3D Contact Engine, Matplotlib (300 DPI) | 3D contacts (HB, SB, $\pi$, Hal, Hyd); co-crystal Tanimoto $T_{\text{IFP}}$; 2D pharmacophore maps |
| **Stage 6** | `src/stage6_md_triage.py` | [OpenMM](https://openmm.org/) 8 (`LangevinMiddleIntegrator`), [GROMACS](https://www.gromacs.org/) 2024 | Three-zone 300 K Langevin MD ($\gamma = 18.0\text{ ps}^{-1}$); `STABLE` / `MODERATE` / `UNSTABLE` triage; Mode C 100-ns export |

## Multi-Target Benchmark Performance (v1.2.0)

AnCore v1.2.0 is validated across five diverse crystallographic benchmark targets from the [RCSB Protein Data Bank](https://www.rcsb.org/) spanning 50 benchmark ligands (25 known high-affinity actives and 25 property-matched ZINC decoys) plus 5 co-crystal reference complexes (`workspaces/`):

| PDB Target | Target Class & Co-Crystal Ligand | Stage 1 Clash Gate ($\min d_{ij}$) | Stage 2 Redocking $\text{RMSD}_{\text{gate}}$ | Stage 4 Static Consensus $\text{ROC-AUC}$ | Stage 6 Post-MD $\text{ROC-AUC}$ | Stage 6 Post-MD Rank #1 Compound |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **[1IEP](https://www.rcsb.org/structure/1IEP)** | BCR-ABL1 Tyrosine Kinase (`STI`, Imatinib) | $2.608\text{ \AA}$ (PASS) | $0.712\text{ \AA}$ (PASS) | $0.680$ | **$0.800$** ($+0.120$) | Imatinib ($\langle\text{RMSD}\rangle = 1.622\text{ \AA}$, `STABLE`) |
| **[3ERT](https://www.rcsb.org/structure/3ERT)** | Estrogen Receptor $\alpha$ (`OHT`, 4-OHT) | $2.267\text{ \AA}$ (PASS) | $1.373\text{ \AA}$ (PASS) | $1.000$ | **$1.000$** ($0.000$) | Raloxifene ($\langle\text{RMSD}\rangle = 1.482\text{ \AA}$, `STABLE`) |
| **[1UYG](https://www.rcsb.org/structure/1UYG)** | HSP90$\alpha$ Chaperone (`PU2`, Purine) | $2.715\text{ \AA}$ (PASS) | $0.761\text{ \AA}$ (PASS) | $0.460$ | **$0.840$** ($+0.380$) | Geldanamycin ($\langle\text{RMSD}\rangle = 1.716\text{ \AA}$, `STABLE`) |
| **[1M17](https://www.rcsb.org/structure/1M17)** | EGFR Tyrosine Kinase (`AQ4`, Erlotinib) | $2.730\text{ \AA}$ (PASS) | $1.234\text{ \AA}$ (PASS) | $0.760$ | **$0.960$** ($+0.200$) | Lapatinib ($\langle\text{RMSD}\rangle = 1.417\text{ \AA}$, `STABLE`) |
| **[1STP](https://www.rcsb.org/structure/1STP)** | Streptavidin (`BTN`, D-Biotin) | $2.445\text{ \AA}$ (PASS) | $0.449\text{ \AA}$ (PASS) | $0.320$ | **$1.000$** ($+0.680$) | D-Biotin ($\langle\text{RMSD}\rangle = 1.448\text{ \AA}$, `STABLE`) |
| **Suite Mean** | **5 Target Classes (50 Ligands)** | **$2.553 \pm 0.194\text{ \AA}$** | **$0.906 \pm 0.388\text{ \AA}$** | **$0.644 \pm 0.264$** | **$0.920 \pm 0.094$** | **100% True Active at Rank #1 (5/5)** |

## Repository Structure

* `src/` — Six-stage Python biophysical pipeline and verification suite
  * `src/stage1_target_prep.py` — Receptor curation, PDBFixer repair, OpenMM AMBER14/GBn2 minimization, steric clash gate, and 5.0 Å pocket geometry
  * `src/stage2_redock_gate.py` — Two-phase co-crystal redocking gate (Phase A local BFGS + Phase B focused Monte Carlo) and symmetry-corrected Kabsch RMSD
  * `src/stage3_ligand_prep.py` — Multi-format ligand ingestion (`.SDF`/`.MOL`/SMILES), RDKit ETKDGv3 3D embedding, MMFF94s relaxation, and ADMET descriptors
  * `src/stage4_screening.py` — Tri-scoring virtual screening (`Vina`, `Vinardo`, `AD4`), Reciprocal Rank Fusion ($\text{RRF}_{k=10}$), and Z-score standardization
  * `src/stage5_interactions.py` — 3D protein-ligand interaction fingerprinting ($T_{\text{IFP}}$) and 300-DPI 2D pharmacophore diagrams
  * `src/stage6_md_triage.py` — Three-zone 300 K OpenMM Langevin relaxation and Mode C 100-ns explicit-solvent OpenMM/GROMACS package builder
  * `src/ancore_custom_runner.py` — Unified CLI/GUI execution runner for custom PDB targets and custom `.SDF` libraries
  * `src/verify_q1_compliance.py` — Automated mathematical and physical gate compliance auditor
* `ancore-tauri/` (and `src-tauri/`, `ui/`) — Native macOS desktop application shell
  * `ancore-tauri/src-tauri/` — Rust backend (`main.rs`, `Cargo.toml`, `tauri.conf.json`) managing asynchronous HTTP asset streaming and subprocess orchestration
  * `ancore-tauri/ui/index.html` — Hardware-accelerated `WKWebView` + [3Dmol.js](https://3dmol.csb.pitt.edu/) 3D molecular viewport and six-stage interactive control panel
  * `ancore-tauri/build_tauri_dmg.sh` — Native macOS `.app` and `.dmg` bundle build script for Apple Silicon (`aarch64-apple-darwin`)
* `workspaces/` — Complete machine-readable benchmark artifacts across `1IEP`, `3ERT`, `1UYG`, `1M17`, and `1STP`
  * Includes per-target minimized receptors, redocked poses, 50-ligand consensus scorecards, 3D IFP contact tables, 300 K Langevin trajectory logs, Mode C 100-ns OpenMM/GROMACS simulation packages, `stage1_suite_summary.json`–`stage6_suite_summary.json`, and `q1_compliance_audit.json`

## Installation & Quick Start

### 1. Conda Environment Setup (macOS Apple Silicon / Linux)

```bash
git clone https://github.com/norseenhosameldeen/AnCore-Scientific-Studio.git
cd AnCore-Scientific-Studio

# Create and activate the Conda environment
conda env create -f environment.yml
conda activate ancore
```

### 2. Running the Six-Stage Pipeline via CLI

```bash
# Run custom target + custom SDF library through all 6 stages
python3 src/ancore_custom_runner.py --pdb 1IEP --workspace workspaces/1IEP

# Verify mathematical and physical gate compliance across all 5 benchmark targets
python3 src/verify_q1_compliance.py
```

### 3. Building & Launching the Native macOS Tauri Studio GUI

```bash
cd ancore-tauri
chmod +x build_tauri_dmg.sh
./build_tauri_dmg.sh
```

## Citation

If you use **AnCore Scientific Studio** in your research, please cite:

> Abdelhakim, N. H., & Senary, A. M. (2026). *AnCore: A Six-Stage Structure-Based Virtual Screening Pipeline Integrating Two-Phase Redocking, Reciprocal Rank Fusion Consensus Scoring, and Three-Zone 300 K Langevin Dynamics Relaxation*. ChemRxiv. https://www.google.com/search?q=https://doi.org/10.26434/chemrxiv.15012010/v1
