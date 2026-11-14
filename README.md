# AnCore Scientific Studio (v1.2.0)

**A Six-Stage Structure-Based Virtual Screening Workbench Integrating Two-Phase Redocking, Reciprocal Rank Fusion Consensus Scoring, and a Contact-Anchored Elastic-Network Pose-Persistence Heuristic**

**Authors:** Norseen H. Abdelhakim and Ahmed M. Senary  
**Affiliation:** Department of Medicinal Chemistry, Faculty of Pharmacy, Galala University, New Galala 43713, Egypt  
**Repository:** [github.com/norseenhosameldeen/AnCore-Scientific-Studio](https://github.com/norseenhosameldeen/AnCore-Scientific-Studio)  
**Current Release:** `v1.2.0` (Six-Stage Pipeline & 5-Target Benchmark Suite)  
**License:** Apache-2.0 License

## Overview

**AnCore Scientific Studio (v1.2.0)** is an open-source macOS desktop workbench and command-line chemoinformatics pipeline for reproducible multi-ligand structure-based virtual screening (SBVS), multi-scoring consensus rescoring, 3D contact fingerprinting, contact-anchored pose-persistence heuristic triage, and turnkey full-atom molecular dynamics (MD) system preparation.

AnCore v1.2.0 decouples a compiled [Rust](https://www.rust-lang.org/) and [Tauri](https://tauri.app/) desktop shell (`ancore-tauri/src-tauri/`) and hardware-accelerated `WKWebView` / [3Dmol.js](https://3dmol.csb.pitt.edu/) viewport (`ancore-tauri/ui/index.html`) from a six-stage Python biophysical engine (`src/`) powered by [PDBFixer](https://github.com/openmm/pdbfixer), [OpenMM](https://openmm.org/), [RDKit](https://www.rdkit.org/), [OpenBabel](https://openbabel.org/), [AutoDock Vina](https://vina.scripps.edu/), and a turnkey Mode C 100-ns explicit-solvent [GROMACS](https://www.gromacs.org/) / OpenMM simulation package exporter. Rather than serializing multi-megabyte coordinate strings over synchronous Python-to-JavaScript IPC pipes, the Rust backend serves workspace coordinate ensembles (`4.19–6.33 MB` of 80-column PDB complexes per target; `71.92 MB` across the 5-target benchmark suite) over an asynchronous localhost HTTP loopback server (`127.0.0.1`) directly into WebGL buffers.

## Version History & Benchmark Reproducibility Note (`v1.0` vs. `v1.2.0`)

All benchmark metrics in this repository (`workspaces/`) and in the `v1.2.0` manuscript match to the third decimal place:

| Benchmark Metric | Initial Prototype (`v1.0` Preprint) | Current Release (`v1.2.0` Workspace & Manuscript) | Methodological Update in `v1.2.0` |
| :--- | :--- | :--- | :--- |
| **Stage 2 Mean Redocking RMSD** | $0.419\text{ \AA}$ (local BFGS only) | **$0.906 \pm 0.388\text{ \AA}$** (Phase B global MC; **$0.524 \pm 0.228\text{ \AA}$** Phase A local BFGS) | Reports unbiased Phase B global Monte Carlo redocking ($\delta_{\text{pad}} = 3.0\text{ \AA}$, `exhaustiveness = 16`) separately from Phase A local crystallographic refinement |
| **Stage 4 Pre-Relaxation Suite ROC-AUC** | $0.728$ | **$0.584$** (`Vina`-only), **$0.648$** ($\text{RRF}_{k=10}$), **$0.680$** ($\overline{\Delta G}_{\text{dock}}$), **$0.696$** (Mean Z-score $\bar{Z}$) | Evaluates [AutoDock Vina](https://vina.scripps.edu/) 1.2 (`Vina`, `Vinardo`, `AD4`) scoring potentials across all 50 torsion-tree ligands with full single-engine ablations |
| **Stage 5 Top-$T_{\text{IFP}}$ Active in `1IEP`** | Nilotinib ($0.789$) | **Imatinib (`ABL_ACT_01`, $T_{\text{IFP}} = 0.833$)**; Nilotinib ($0.610$); Ponatinib ($0.516$) | Computed directly from 3D heavy-atom radial contact sets against the minimized `1IEP` (`STI_REF`) crystal complex |
| **Stage 6 Heuristic & Control ROC-AUC** | $1.000$ (`5.0 ns` label) | **$0.864$** (Co-Crystal-Blind Contact+Dock); **$0.912$** (LOTO CV); **$0.920 \pm 0.094$** ($S_{\text{relax}}$, 1,500 steps / $3.0\text{ ps}$) | Renames Mode B to a contact-anchored elastic-network heuristic ($1,500\text{ steps} = 3.0\text{ ps}$ at $\Delta t = 2.0\text{ fs}$), adds circularity/LOTO controls, and reports unrestrained full-atom AMBER14/GBn2+GAFF minimization (mean RMSD $= 0.440 \pm 0.131\text{ \AA}$) |
| **Stage 6 Discrete Cutoff Verdicts** | 25/25 decoys rejected | **Decoys:** `22/25 UNSTABLE`, `3/25 METASTABLE`, `0/25 STABLE`; **Actives:** `6/25 STABLE`, `9/25 METASTABLE`, `10/25 UNSTABLE` | Discloses the 40% false-negative rate of hard cutoff gating on multi-rotor actives (`Imatinib`, `4-OHT`, `Erlotinib`), motivating Mode C full-atom explicit-solvent simulation |

## Six-Stage Pipeline Architecture

| Stage | Module Script | Core Scientific Engines | Quality-Control Gate & Output Metrics |
| :--- | :--- | :--- | :--- |
| **Stage 1** | `src/stage1_target_prep.py` | [PDBFixer](https://github.com/openmm/pdbfixer), [OpenMM](https://openmm.org/) (`amber14-all` + `gbn2`), [OpenBabel](https://openbabel.org/) | pH 7.4 protonation; 200-step L-BFGS minimization; steric clash gate $\min d_{ij} \ge 2.10\text{ \AA}$; 5.0 Å pocket box |
| **Stage 2** | `src/stage2_redock_gate.py` | [AutoDock Vina](https://vina.scripps.edu/) 1.2, SciPy linear-sum assignment RMSD | Phase A local BFGS + Phase B focused MC ($\delta_{\text{pad}} = 3.0\text{ \AA}$, `exhaustiveness = 16`, seed = 42); symmetry-corrected $\text{RMSD} \le 2.00\text{ \AA}$ |
| **Stage 3** | `src/stage3_ligand_prep.py` | [RDKit](https://www.rdkit.org/) (`ETKDGv3`, `MMFF94s`), [OpenBabel](https://openbabel.org/) | `.SDF`/`.MOL`/SMILES import; 500-iter MMFF94s relaxation; MW, cLogP, TPSA, HBD, HBA, $N_{\text{rot}}$, QED |
| **Stage 4** | `src/stage4_screening.py` | [AutoDock Vina](https://vina.scripps.edu/) (`Vina`, `Vinardo`, `AD4` rescoring) | Discovery box ($\delta_{\text{pad}} = 4.0\text{ \AA}$, `exhaustiveness = 8`); $\text{RRF}_{k=10}$; mean Z-score $\bar{Z}$; single-scorer ablations |
| **Stage 5** | `src/stage5_interactions.py` | AnCore 3D Radial Contact Profiler, Matplotlib (300 DPI) | 3D radial contacts (HB $\le 3.5\text{ \AA}$, SB $\le 4.2\text{ \AA}$, $\pi \le 4.3\text{ \AA}$, Hal $\le 3.9\text{ \AA}$, Hyd $\le 4.0\text{ \AA}$); $T_{\text{IFP}}$; 2D maps |
| **Stage 6** | `src/stage6_md_triage.py` | [OpenMM](https://openmm.org/) 8 (`LangevinMiddleIntegrator`), [GROMACS](https://www.gromacs.org/) 2024 | Contact-anchored elastic-network heuristic (Mode B: 1,500 steps / $3.0\text{ ps}$, $\gamma = 18.0\text{ ps}^{-1}$, $\Delta t = 2.0\text{ fs}$); full-atom AMBER14/GBn2+GAFF baseline & Mode C 100-ns export |

## Multi-Target Benchmark Performance & Circularity Controls (`workspaces/`)

| PDB Target | Co-Crystal Ref. | Stage 2 Redock RMSD (Phase B / A) | Stage 4 Static AUC (`Vina` / `RRF` / $\bar{Z}$) | Co-Crystal-Blind Contact+Dock AUC | LOTO Cross-Validated AUC | Stage 6 Heuristic $S_{\text{relax}}$ AUC | Top Lead Full-Atom AMBER14/GBn2+GAFF RMSD ($E_{\min}$) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **[1IEP](https://www.rcsb.org/structure/1IEP)** | `STI` (Imatinib) | $0.712\text{ \AA}$ / $0.337\text{ \AA}$ | $0.680$ / $0.680$ / $0.680$ | $0.840$ | $0.840$ | **$0.800$** | Nilotinib: **$0.585\text{ \AA}$** ($-9,597.8\text{ kcal/mol}$) |
| **[3ERT](https://www.rcsb.org/structure/3ERT)** | `OHT` (4-OHT) | $1.373\text{ \AA}$ / $0.488\text{ \AA}$ | $1.000$ / $1.000$ / $1.000$ | $1.000$ | $0.960$ | **$1.000$** | Raloxifene: **$0.405\text{ \AA}$** ($-7,958.0\text{ kcal/mol}$) |
| **[1UYG](https://www.rcsb.org/structure/1UYG)** | `PU2` (Purine) | $0.761\text{ \AA}$ / $0.640\text{ \AA}$ | $0.440$ / $0.480$ / $0.600$ | $0.800$ | $0.800$ | **$0.840$** | Luminespib: **$0.294\text{ \AA}$** ($-6,958.1\text{ kcal/mol}$) |
| **[1M17](https://www.rcsb.org/structure/1M17)** | `AQ4` (Erlotinib) | $1.234\text{ \AA}$ / $0.855\text{ \AA}$ | $0.600$ / $0.760$ / $0.800$ | $0.800$ | $0.960$ | **$0.960$** | Lapatinib: **$0.569\text{ \AA}$** ($-15,455.9\text{ kcal/mol}$) |
| **[1STP](https://www.rcsb.org/structure/1STP)** | `BTN` (D-Biotin) | $0.449\text{ \AA}$ / $0.302\text{ \AA}$ | $0.200$ / $0.320$ / $0.400$ | $0.880$ | $1.000$ | **$1.000$** | D-Biotin: **$0.348\text{ \AA}$** ($-4,035.2\text{ kcal/mol}$) |
| **Suite Mean** | **5 Targets (50 Ligands)** | **$0.906\text{ \AA}$ / $0.524\text{ \AA}$** | **$0.584$ / $0.648$ / $0.696$** | **$0.864 \pm 0.083$** | **$0.912 \pm 0.087$** | **$0.920 \pm 0.094$** | **Mean: $0.440 \pm 0.131\text{ \AA}$ (5/5 $\le 0.60\text{ \AA}$)** |

* **Pooled 50-Compound 95% Bootstrap Confidence Intervals (2,000 Resamples):** `Vina-only: 0.586 [0.406, 0.750]`, `RRF: 0.681 [0.522, 0.826]`, `Mean-Z: 0.694 [0.537, 0.837]`, `Mode B Heuristic <RMSD>: 0.872 [0.767, 0.955]`, `Mode B Heuristic S_relax: 0.895 [0.798, 0.967]`.
* **Shuffled-Contact Negative Control (200 Random Permutations/Target):** Collapses suite mean ROC-AUC from $0.920$ to **$0.594 \pm 0.136$** (`0.676`, `0.792`, `0.510`, `0.534`, `0.460`).
* **Unbiased Physical AMBER/GAFF LJ+Polar Relaxation Control (Zero `I_anchor`, Zero Heuristic Contact Springs):** Achieves a suite mean ROC-AUC of **$0.688 \pm 0.219$** (`+0.104` over `Vina-only`).

## Installation & Quick Start

```bash
git clone https://github.com/norseenhosameldeen/AnCore-Scientific-Studio.git
cd AnCore-Scientific-Studio
conda env create -f environment.yml
conda activate ancore
python3 src/ancore_custom_runner.py --pdb 1STP --workspace workspaces/1STP
```

## Citation

> Abdelhakim, N. H., & Senary, A. M. (2026). *AnCore: A Six-Stage Structure-Based Virtual Screening Pipeline Integrating Two-Phase Redocking, Reciprocal Rank Fusion Consensus Scoring, and a Contact-Anchored Elastic-Network Pose-Persistence Heuristic*. ChemRxiv. https://www.google.com/search?q=https://doi.org/10.26434/chemrxiv.15012010/v1
