# AnCore Scientific Studio (v1.2.0)

**A Six-Stage Structure-Based Virtual Screening Workbench Integrating Two-Phase Redocking, Reciprocal Rank Fusion Consensus Rescoring, and Three-Zone 300 K Langevin Dynamics Pose Triage**

**Authors:** Norseen H. Abdelhakim and Ahmed M. Senary  
**Affiliation:** Department of Medicinal Chemistry, Faculty of Pharmacy, Galala University, New Galala 43713, Egypt  
**Repository:** [github.com/norseenhosameldeen/AnCore-Scientific-Studio](https://github.com/norseenhosameldeen/AnCore-Scientific-Studio)  
**Preprint DOI:** [10.26434/chemrxiv.15012010/v1](https://www.google.com/search?q=https://doi.org/10.26434/chemrxiv.15012010/v1)  
**Current Release:** `v1.2.0` (Full Six-Stage Physics & 5-Target Benchmark Suite)  
**License:** MIT License

## Overview

**AnCore Scientific Studio (v1.2.0)** is an open-source macOS desktop workbench and command-line chemoinformatics pipeline for reproducible multi-ligand structure-based virtual screening (SBVS) and post-docking pose stability triage.

AnCore v1.2.0 integrates a compiled [Rust](https://www.rust-lang.org/) and [Tauri](https://tauri.app/) desktop shell (`ancore-tauri/src-tauri/`) and hardware-accelerated `WKWebView` / [3Dmol.js](https://3dmol.csb.pitt.edu/) viewport (`ancore-tauri/ui/index.html`) with a six-stage Python biophysical engine (`src/`) powered by [PDBFixer](https://github.com/openmm/pdbfixer), [OpenMM](https://openmm.org/), [RDKit](https://www.rdkit.org/), [OpenBabel](https://openbabel.org/), [AutoDock Vina](https://vina.scripps.edu/), and a turnkey Mode C 100-ns explicit-solvent [GROMACS](https://www.gromacs.org/) / OpenMM simulation package exporter. Rather than serializing multi-megabyte coordinate strings over synchronous Python-to-JavaScript IPC pipes, the Rust backend serves workspace coordinate ensembles (`4.19–6.33 MB` of 80-column PDB complexes per target; `71.92 MB` across the 5-target benchmark suite) over an asynchronous localhost HTTP loopback server directly into WebGL buffers.

## Version History & Benchmark Reproducibility Note (`v1.0` vs. `v1.2.0`)

To ensure complete transparency between the initial ChemRxiv preprint (`v1.0`) and the revised full-length manuscript (`v1.2.0`), all benchmark metrics in this repository and in the manuscript are synchronized directly to the machine-readable JSON/CSV manifests in `workspaces/`:

| Benchmark Metric | Initial Prototype (`v1.0` Preprint) | Current Release (`v1.2.0` Workspace & Manuscript) | Methodological Reason for Update in `v1.2.0` |
| :--- | :--- | :--- | :--- |
| **Stage 2 Mean Redocking RMSD** | $0.419\text{ \AA}$ (local BFGS only) | **$0.906 \pm 0.388\text{ \AA}$** (Phase B global MC; **$0.524 \pm 0.228\text{ \AA}$** Phase A local BFGS) | Reports unbiased Phase B global Monte Carlo redocking ($\delta_{\text{pad}} = 3.0\text{ \AA}$, `exhaustiveness = 16`) separately from Phase A local crystallographic refinement |
| **Stage 4 Pre-Relaxation Suite ROC-AUC** | $0.728$ | **$0.584$** (`Vina`-only), **$0.648$** ($\text{RRF}_{k=10}$), **$0.644$** ($\overline{\Delta G}_{\text{dock}}$), **$0.696$** (Mean Z-score $\bar{Z}$) | Evaluates real [AutoDock Vina](https://vina.scripps.edu/) 1.2 (`Vina`, `Vinardo`, `AD4`) scoring potentials across all 50 torsion-tree ligands |
| **Stage 5 Top-$T_{\text{IFP}}$ Active in `1IEP`** | Nilotinib ($0.789$) | **Imatinib (`ABL_ACT_01`, $T_{\text{IFP}} = 0.833$)**; Nilotinib ($0.610$); Ponatinib ($0.516$) | Computed directly from 3D heavy-atom contact sets against the minimized `1IEP` (`STI_REF`) crystal complex |
| **Stage 6 Post-Relaxation Suite ROC-AUC** | $1.000$ (5/5 targets) | **$0.896 \pm 0.096$** (raw unweighted $\langle\text{RMSD}\rangle$); **$0.920 \pm 0.094$** ($S_{\text{relax}}$: range $0.800\text{--}1.000$) | Replaces heuristic interpolation with 1,500-step ($3.0\text{ ps}$) 300 K OpenMM Langevin trajectories (`stage6_md_triage.py`) |
| **Stage 6 Decoy Triage Classification** | 25/25 rejected | **22/25 (`88.0%`) `UNSTABLE (REJECT)`**, **3/25 (`12.0%`) `METASTABLE (WARN)`**, **0/25 (`0.0%`) `STABLE (PASS)`** | Matches Appendix Tables S1–S5 exactly (`ZINC_Decoy_A2`, `ZINC_Decoy_S4`, and `ZINC_Decoy_S5` are `METASTABLE`) |

## Six-Stage Pipeline Architecture

| Stage | Module Script | Core Scientific Engines | Quality-Control Gate & Output Metrics |
| :--- | :--- | :--- | :--- |
| **Stage 1** | `src/stage1_target_prep.py` | [PDBFixer](https://github.com/openmm/pdbfixer), [OpenMM](https://openmm.org/) (`amber14-all` + `gbn2`), [OpenBabel](https://openbabel.org/) | pH 7.4 protonation; 200-step L-BFGS minimization; steric clash gate $\min d_{ij} \ge 2.10\text{ \AA}$; 5.0 Å pocket box |
| **Stage 2** | `src/stage2_redock_gate.py` | [AutoDock Vina](https://vina.scripps.edu/) 1.2, SciPy linear-sum assignment RMSD | Phase A local BFGS + Phase B focused MC ($\delta_{\text{pad}} = 3.0\text{ \AA}$, `exhaustiveness = 16`); gate $\text{RMSD} \le 2.00\text{ \AA}$ |
| **Stage 3** | `src/stage3_ligand_prep.py` | [RDKit](https://www.rdkit.org/) (`ETKDGv3`, `MMFF94s`), [OpenBabel](https://openbabel.org/) | `.SDF`/`.MOL`/SMILES import; 500-iter MMFF94s relaxation; MW, cLogP, TPSA, HBD, HBA, $N_{\text{rot}}$, QED |
| **Stage 4** | `src/stage4_screening.py` | [AutoDock Vina](https://vina.scripps.edu/) (`Vina`, `Vinardo`, `AD4` rescoring) | Discovery box ($\delta_{\text{pad}} = 4.0\text{ \AA}$, `exhaustiveness = 8`); $\text{RRF}_{k=10}$; mean Z-score $\bar{Z}$; single-scorer ablations |
| **Stage 5** | `src/stage5_interactions.py` | AnCore 3D Radial Contact Profiler, Matplotlib (300 DPI) | 3D radial contacts (HB $\le 3.5\text{ \AA}$, SB $\le 4.2\text{ \AA}$, $\pi \le 4.3\text{ \AA}$, Hal $\le 3.9\text{ \AA}$, Hyd $\le 4.0\text{ \AA}$); $T_{\text{IFP}}$; 2D maps |
| **Stage 6** | `src/stage6_md_triage.py` | [OpenMM](https://openmm.org/) 8 (`LangevinMiddleIntegrator`), [GROMACS](https://www.gromacs.org/) 2024 | Three-zone 300 K Langevin relaxation (1,500 steps, $\gamma = 18.0\text{ ps}^{-1}$, $\Delta t = 2.0\text{ fs}$); stability score $S_{\text{relax}}$; Mode C 100-ns export |

## Multi-Target Benchmark Performance (`workspaces/`)

Validated across five crystallographic targets from the [RCSB Protein Data Bank](https://www.rcsb.org/) spanning 50 benchmark ligands (25 verified actives and 25 property-matched ZINC decoys) plus 5 co-crystal reference complexes:

| PDB Target | Target Class & Co-Crystal Ref. | Stage 1 $\min d_{ij}$ | Stage 2 Redock RMSD (Phase B / Phase A) | Stage 4 Static AUC (`Vina` / `RRF` / $\bar{Z}$) | Stage 6 Post-MD AUC (Raw $\langle\text{RMSD}\rangle$ / $S_{\text{relax}}$) | Stage 6 Post-MD Rank #1 Lead ($\langle\text{RMSD}\rangle$, $\text{HB}_{\text{occ}}$, Verdict) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **[1IEP](https://www.rcsb.org/structure/1IEP)** | BCR-ABL1 Kinase (`STI`, Imatinib) | $2.608\text{ \AA}$ | $0.712\text{ \AA}$ / $0.337\text{ \AA}$ | $0.680$ / $0.680$ / $0.680$ | $0.760$ / **$0.800$** | Nilotinib (`ABL_ACT_02`, $2.562\text{ \AA}$, $0.0\%$, `METASTABLE`) |
| **[3ERT](https://www.rcsb.org/structure/3ERT)** | Estrogen Receptor $\alpha$ (`OHT`, 4-OHT) | $2.267\text{ \AA}$ | $1.373\text{ \AA}$ / $0.488\text{ \AA}$ | $1.000$ / $1.000$ / $1.000$ | $0.920$ / **$1.000$** | Raloxifene (`ERA_ACT_02`, $2.134\text{ \AA}$, $20.0\%$, `METASTABLE`) |
| **[1UYG](https://www.rcsb.org/structure/1UYG)** | HSP90$\alpha$ Chaperone (`PU2`, Purine) | $2.715\text{ \AA}$ | $0.761\text{ \AA}$ / $0.640\text{ \AA}$ | $0.440$ / $0.480$ / $0.600$ | $0.840$ / **$0.840$** | Luminespib (`HSP_ACT_04`, $1.364\text{ \AA}$, $80.5\%$, `STABLE`) |
| **[1M17](https://www.rcsb.org/structure/1M17)** | EGFR Kinase (`AQ4`, Erlotinib) | $2.730\text{ \AA}$ | $1.234\text{ \AA}$ / $0.855\text{ \AA}$ | $0.600$ / $0.760$ / $0.800$ | $0.960$ / **$0.960$** | Lapatinib (`EGF_ACT_04`, $1.257\text{ \AA}$, $74.4\%$, `STABLE`) |
| **[1STP](https://www.rcsb.org/structure/1STP)** | Streptavidin (`BTN`, D-Biotin) | $2.445\text{ \AA}$ | $0.449\text{ \AA}$ / $0.302\text{ \AA}$ | $0.200$ / $0.320$ / $0.400$ | $1.000$ / **$1.000$** | D-Biotin (`STP_ACT_01`, $1.795\text{ \AA}$, $83.8\%$, `STABLE`) |
| **Suite Mean** | **5 Target Classes (50 Ligands)** | **$2.553\text{ \AA}$** | **$0.906\text{ \AA}$ / $0.524\text{ \AA}$** | **$0.584$ / $0.648$ / $0.696$** | **$0.896$ / $0.920$** | **Active at Rank #1 on 5/5 Targets (3 `STABLE`, 2 `METASTABLE`)** |

### Summary of Stage 6 Stability Verdicts Across All 50 Screened Complexes

* **Verified Active Inhibitors ($N = 25$):** `6/25 (24.0%) STABLE (PASS)`, `9/25 (36.0%) METASTABLE (WARN)`, and `10/25 (40.0%) UNSTABLE (REJECT)`. Notably, in deep or water-bridged pockets (`1IEP`, `3ERT`) and when upstream Stage 4 discovery-box sampling (`exhaustiveness = 8`, $\delta_{\text{pad}} = 4.0\text{ \AA}$) mis-docks high-rotor inhibitors such as Erlotinib (`EGF_ACT_01`, $N_{\text{rot}} = 10$, Stage 4 Rank #10), unsolvated coarse-grained relaxation penalizes flexible solvent-exposed tails (`Imatinib`: $3.580\text{ \AA}$; `4-Hydroxytamoxifen`: $4.167\text{ \AA}$; `Erlotinib`: $4.818\text{ \AA}$), motivating the automated **Mode C** 100-ns explicit-solvent export.
* **Property-Matched ZINC Decoys ($N = 25$):** `0/25 (0.0%) STABLE (PASS)` ($100\%$ precision for the `STABLE` tier), `3/25 (12.0%) METASTABLE (WARN)` (`ZINC_Decoy_A2` in `1IEP`, `ZINC_Decoy_S4` and `ZINC_Decoy_S5` in `1STP`), and `22/25 (88.0%) UNSTABLE (REJECT)`.

## Repository Structure

* `src/` — Six-stage Python biophysical pipeline and verification suite
  * `src/stage1_target_prep.py` — Receptor curation, PDBFixer repair, OpenMM AMBER14/GBn2 minimization, steric clash gate, and 5.0 Å pocket geometry
  * `src/stage2_redock_gate.py` — Two-phase co-crystal redocking gate (Phase A local BFGS + Phase B focused Monte Carlo) and linear-sum assignment RMSD
  * `src/stage3_ligand_prep.py` — Multi-format ligand ingestion (`.SDF`/`.MOL`/SMILES), RDKit ETKDGv3 3D embedding, MMFF94s relaxation, and ADMET descriptors
  * `src/stage4_screening.py` — Tri-scoring virtual screening (`Vina`, `Vinardo`, `AD4`), Reciprocal Rank Fusion ($\text{RRF}_{k=10}$), and Z-score standardization
  * `src/stage5_interactions.py` — 3D radial protein-ligand interaction fingerprinting ($T_{\text{IFP}}$) and 300-DPI 2D pharmacophore diagrams
  * `src/stage6_md_triage.py` — Three-zone 300 K OpenMM Langevin relaxation (1,500 steps, $\Delta t = 2.0\text{ fs}$) and Mode C 100-ns explicit-solvent OpenMM/GROMACS package builder
  * `src/ancore_custom_runner.py` — Unified CLI/GUI execution runner for custom PDB targets and custom `.SDF` libraries
  * `src/verify_q1_compliance.py` — Automated mathematical and physical gate compliance auditor
* `ancore-tauri/` (and `src-tauri/`, `ui/`) — Native macOS desktop application shell
  * `ancore-tauri/src-tauri/` — Rust backend (`main.rs`, `Cargo.toml`, `tauri.conf.json`) managing asynchronous localhost HTTP asset streaming and subprocess orchestration
  * `ancore-tauri/ui/index.html` — Hardware-accelerated `WKWebView` + [3Dmol.js](https://3dmol.csb.pitt.edu/) 3D molecular viewport and six-stage interactive control panel
  * `ancore-tauri/build_tauri_dmg.sh` — Native macOS `.app` and `.dmg` bundle build script for Apple Silicon (`aarch64-apple-darwin`)
* `workspaces/` — Complete machine-readable benchmark artifacts across `1IEP`, `3ERT`, `1UYG`, `1M17`, and `1STP` (`stage1_suite_summary.json`–`stage6_suite_summary.json`, `q1_compliance_audit.json`, pre-merged 80-column PDB complexes, trajectory CSVs, and Mode C 100-ns OpenMM/GROMACS packages)

## Installation & Usage

```bash
git clone https://github.com/norseenhosameldeen/AnCore-Scientific-Studio.git
cd AnCore-Scientific-Studio

# 1. Create and activate the Conda environment
conda env create -f environment.yml
conda activate ancore

# 2. Verify all 5 benchmark targets and audit mathematical/gate compliance
python3 src/verify_q1_compliance.py

# 3. Run the 6-stage pipeline on a target workspace
python3 src/ancore_custom_runner.py --pdb 1STP --workspace workspaces/1STP
```

## Citation

> Abdelhakim, N. H., & Senary, A. M. (2026). *AnCore: A Six-Stage Structure-Based Virtual Screening Pipeline Integrating Two-Phase Redocking, Reciprocal Rank Fusion Consensus Scoring, and Three-Zone 300 K Langevin Dynamics Relaxation*. ChemRxiv. https://www.google.com/search?q=https://doi.org/10.26434/chemrxiv.15012010/v1
