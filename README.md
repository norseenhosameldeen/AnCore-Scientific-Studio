# AnCore Scientific Studio (v1.2.0)

**A six-stage structure-based virtual screening workbench integrating two-phase redocking, Reciprocal Rank Fusion consensus scoring, and a contact-anchored elastic-network pose-persistence heuristic**

**Authors:** Norseen H. Abdelhakim and Ahmed M. Senary
**Affiliation:** Department of Medicinal Chemistry, Faculty of Pharmacy, Galala University, New Galala 43713, Egypt
**Release:** `v1.2.0` | **License:** Apache-2.0 | **Preprint:** [ChemRxiv 10.26434/chemrxiv.15012010](https://doi.org/10.26434/chemrxiv.15012010)
**Archived snapshot (Zenodo):** _DOI to be added_

---

## 1. What this is (and is not)

AnCore is an open-source research workbench with two parts:

* a **macOS desktop shell** (Rust/Tauri 1.6 + `WKWebView` + 3Dmol.js) in `ancore-tauri/` that loads coordinates and stage tables from a local workspace over a `127.0.0.1` HTTP loopback server;
* a **six-stage Python engine** in `src/` (PDBFixer, OpenMM 8, RDKit, OpenBabel/Meeko, AutoDock Vina 1.2) that writes JSON manifests and CSV audit tables for every stage.

**Scope and status.** The paper's results were produced on macOS (Apple Silicon, 8-core M1). The desktop application is macOS-only. The Python engine depends on the conda environment below and has not been tested on other platforms by the authors.

**Stage 6 Mode B is not molecular dynamics.** It is a short (1,500 steps, 3.0 ps) contact-anchored elastic-network heuristic whose restraints are scaled by the same Stage 5 contact counts it is later used to rank. It is a ranking heuristic, not a free-energy method and not a validated pose-stability predictor. Mode C only *exports* 100-ns explicit-solvent packages; the 100-ns simulations themselves have not been run or validated here.

---

## 2. Version history and benchmark reproducibility note (v1.0 → v1.2.0)

The initial prototype (v1.0, ChemRxiv v1) reported different numbers, including a Stage 6 ROC-AUC of 1.000 and a "5.0 ns" label for what is a 3.0 ps heuristic. Those results are **superseded**. v1.2.0 re-implemented Stages 2, 4 and 6, renamed Mode B, and added controls. All v1.2.0 metrics in `workspaces/` correspond to the v1.2.0 manuscript.

| Benchmark metric | v1.0 (superseded) | v1.2.0 (this repository) | What changed |
| :--- | :--- | :--- | :--- |
| Stage 2 mean redocking RMSD | 0.419 Å (local BFGS only) | **0.906 ± 0.388 Å** (Phase B global MC); **0.524 ± 0.228 Å** (Phase A local BFGS) | Phase B (unbiased global search) reported separately from Phase A (crystal-seeded refinement) |
| Stage 4 suite ROC-AUC | 0.728 | **0.584** Vina · **0.648** RRF (k = 10) · **0.680** mean ΔG · **0.696** mean Z | Tri-scoring with single-engine ablations |
| Stage 5 top-T_IFP active, 1IEP | Nilotinib (0.789) | **Imatinib (ABL_ACT_01), T_IFP = 0.833** | T_IFP computed from 3D heavy-atom contacts against the minimized STI reference |
| Stage 6 ROC-AUC | 1.000 ("5.0 ns") | **0.864** co-crystal-blind contact + dock · **0.912** LOTO · **0.920 ± 0.094** S_relax | Mode B relabelled as a 3.0 ps heuristic; circularity controls added |
| Stage 6 discrete verdicts | 25/25 decoys rejected | Decoys: 22 UNSTABLE / 3 METASTABLE / 0 STABLE. Actives: 6 STABLE / 9 METASTABLE / **10 UNSTABLE** | The 40% false-negative rate on actives is disclosed |

> **Note on cached or older copies.** Older README text, figures, or the ChemRxiv v1 preprint may show the v1.0 numbers above. The tag `v1.2.0` and the files under `workspaces/` are authoritative.

---

## 3. Pipeline

| Stage | Script | Engines | Gate / key parameters |
| :-- | :-- | :-- | :-- |
| 1 | `src/stage1_target_prep.py` | PDBFixer, OpenMM (`amber14-all` + `implicit/gbn2`), OpenBabel | pH 7.4 protonation; 200-iteration L-BFGS; steric clash gate min d_ij ≥ 2.10 Å; 5.0 Å pocket mapping |
| 2 | `src/stage2_redock_gate.py` | AutoDock Vina 1.2, SciPy `linear_sum_assignment` | Phase A local BFGS + Phase B focused MC (δ_pad = 3.0 Å, exhaustiveness 16, seed 42); symmetry-corrected RMSD ≤ 2.00 Å |
| 3 | `src/stage3_ligand_prep.py` | RDKit (ETKDGv3, MMFF94s), OpenBabel/Meeko | `.sdf`/`.mol`/SMILES; 500-iteration MMFF94s; MW, cLogP, TPSA, HBD, HBA, N_rot, QED |
| 4 | `src/stage4_screening.py` | Vina, Vinardo, AutoDock4 | Discovery box (δ_pad = 4.0 Å, exhaustiveness 8, seed 42); RRF (k = 10); mean Z-score |
| 5 | `src/stage5_interactions.py` | AnCore 3D radial contact engine | HB ≤ 3.5, SB ≤ 4.2, π ≤ 4.3, halogen ≤ 3.9, hydrophobic ≤ 4.0 Å; T_IFP vs. co-crystal reference. Heavy-atom radial cutoffs only; no donor–H–acceptor angle filtering |
| 6 | `src/stage6_md_triage.py` | OpenMM 8 (`LangevinMiddleIntegrator`), GROMACS 2024 exporter | Mode A: 100-step unrestrained full-atom minimization. Mode B: 1,500 steps (50 frames × 30 steps, Δt = 2.0 fs, 3.0 ps), γ = 18.0 ps⁻¹, 300 K, seed 42. Mode C: 100-ns package export |

**Mode B implementation details (as coded in v1.2.0).**

* Receptor atoms are partitioned into three radial zones about the pocket centre: inner (≤ 8.0 Å, harmonic restraints k = 1,450 kJ mol⁻¹ nm⁻² on backbone N/CA/C/O and 680 on side-chain heavy atoms), buffer (8–12 Å, k = 2,092), and massless outer scaffold (> 12 Å). The receptor is truncated at **15.0 Å** from the pocket centre at the call site in `stage6_md_triage.py` (the function default of 16.0 Å is overridden).
* Ligand heavy atoms form an elastic network; Stage 5 polar and hydrophobic contacts become harmonic springs; a continuous outward pull `f_desolv = max(14, 290 − 15.5·I_anchor)` and confinement `k_conf = 180 + 18.5·I_anchor` act on the ligand.
* `I_anchor` is built from Stage 5 contact counts (`N_polar`, `N_ref`, `T_IFP`), ligand HBD/HBA, heavy-atom count and a **rigid-steroid term** that applies when N_rot = 0 and E_Vina ≤ −9.0 kcal/mol (+3.5 in `I_anchor`, and +2.5 in the effective-anchor count of S_relax). In the benchmark this term affects only 17β-estradiol (3ERT).
* **Parameter provenance.** All Mode B constants were set during development while inspecting results on the five benchmark targets. They are not independent of the benchmark (see Section 5).
* **Time axis.** `time_ns_series` in the exported trajectory CSVs/figures is a legacy plotting scale (5.0 × frame / n_frames). The integration itself is 3.0 ps.

---

## 4. Benchmark (development set)

**Composition.** 5 targets × 10 ligands = 50 screened complexes (25 literature actives, 25 presumed-inactive decoys), plus 5 co-crystal references (55 complexes). Actives are hand-selected literature compounds (potencies in manuscript Table S6). Decoys are generic ZINC-style compounds that are *presumed* inactive, **not** DUD-E-generated, and **not strictly property-matched** (mean N_rot 6.32 vs. 4.08, actives vs. decoys). The five co-crystal ligands (or their close analogs) are also among the actives, and the co-crystal reference defines `N_ref` and `T_IFP`.

| PDB | Co-crystal ref. | Stage 2 redock RMSD (Phase B / A) | Stage 4 AUC (Vina / RRF / Z̄) | Co-crystal-blind contact + dock AUC | LOTO AUC | Stage 6 S_relax AUC | Mode A top-lead RMSD (E_min) |
| :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| 1IEP | STI (imatinib) | 0.712 / 0.337 Å | 0.680 / 0.680 / 0.680 | 0.840 | 0.840 | 0.800 | Nilotinib 0.585 Å (−9,597.8 kcal/mol) |
| 3ERT | OHT (4-OHT) | 1.373 / 0.488 Å | 1.000 / 1.000 / 1.000 | 1.000 | 0.960 | 1.000 | Raloxifene 0.405 Å (−7,958.0) |
| 1UYG | PU2 (purine) | 0.761 / 0.640 Å | 0.440 / 0.480 / 0.600 | 0.800 | 0.800 | 0.840 | Luminespib 0.294 Å (−6,958.1) |
| 1M17 | AQ4 (erlotinib) | 1.234 / 0.855 Å | 0.600 / 0.760 / 0.800 | 0.800 | 0.960 | 0.960 | Lapatinib 0.569 Å (−15,455.9) |
| 1STP | BTN (D-biotin) | 0.449 / 0.302 Å | 0.200 / 0.320 / 0.400 | 0.880 | 1.000 | 1.000 | D-Biotin 0.348 Å (−4,035.2) |
| **Mean** | 5 targets | **0.906 / 0.524 Å** | **0.584 / 0.648 / 0.696** | **0.864 ± 0.083** | **0.912 ± 0.087** | **0.920 ± 0.094** | **0.440 ± 0.131 Å** |

* Within-target stratified bootstrap (2,000 resamples, 5 actives + 5 decoys resampled per target): Vina-only 0.584 [0.440, 0.736]; mean Z 0.696 [0.552, 0.840]; S_relax 0.920 [0.836, 0.992]. Reported p-values are one-sided paired bootstrap p-values versus Vina-only.
* Mode B verdict audit: decoys 22 UNSTABLE / 3 METASTABLE / 0 STABLE; actives 6 STABLE / 9 METASTABLE / 10 UNSTABLE (**40% false-negative rate**). Do **not** use the Mode B UNSTABLE label as a rejection filter for flexible multi-rotor ligands.
* The shuffled control in `run_ablation_and_controls.py` permutes the Mode B ⟨RMSD⟩ values across ligands within each target (200 permutations, seed 42), then rescoring. It tests whether the pairing of drift with ligand identity matters. It does not re-run Mode B with permuted contacts.

---

## 5. Limitations and honest interpretation

1. **Small, non-standard benchmark.** AUC rests on 5 actives vs. 5 decoys per target; EF20% takes only the values {0, 1, 2}. DUD-E/LIT-PCBA-scale evaluation (≈1:50 ratios) has not been done.
2. **In-sample Stage 6 results.** Mode B constants were set while viewing these five targets. LOTO refits only linear weights on features designed using all five targets, so it does not remove this dependence. Treat the 0.920 (and 0.864/0.912) figures as in-sample, descriptive results.
3. **Input/output overlap in Mode B.** Restraint strengths depend on Stage 5 contact counts; ranking by Mode B outputs partly re-measures those counts. Most of the gain over Vina-only is already reached by the co-crystal-blind contact-count control (0.864).
4. **No ligand-property or other trivial baselines** have been reported yet (e.g., HBD/HBA/N_rot-only scoring).
5. **Single docking seed** (42) for all Vina runs. No replicate-seed confidence intervals.
6. **Mode A** (100-step minimization) shows only local relaxation of docked poses; it is not evidence of pose correctness.
7. **Stage 2 gate.** Phase A starts from the crystal pose and is pooled with Phase B poses, so a pass is partly seeded by the answer. The 1M17 naive-box failure (5.549 Å) is a single-seed observation.
8. **Stage 5** uses heavy-atom radial cutoffs without angle filtering.

---

## 6. Reproducibility status (v1.2.0)

| Item | Location | Status |
| :-- | :-- | :-- |
| Per-target stage manifests and CSV/PDB outputs | `workspaces/<PDB>/` | Included |
| Suite summaries (Stages 1–6) | `workspaces/stage{1..6}_suite_summary.json` | Included; Stage 2, 4 and 6 values cross-checked against the manuscript |
| Single-engine/RRF ablations, N_polar and T_IFP controls, co-crystal-blind control, shuffled-drift control, LOTO, stratified bootstrap | `run_ablation_and_controls.py` (repository root) | Included; run `python3 run_ablation_and_controls.py --workspaces-dir workspaces` |
| Integrity checks on workspaces/source | `src/verify_q1_compliance.py` → `workspaces/q1_compliance_audit.json` | Included |
| 1M17 focused-grid (δ_pad = 3.0 Å, exhaustiveness 16) Stage 4–6 rerun | — | **Artifacts not yet in `workspaces/`**; to be added |
| Unbiased AMBER/GAFF 12–6/10–12 relaxation control (AUC 0.688) | — | **Script/outputs not yet located**; to be added |
| Pooled 50-compound bootstrap CIs (e.g., S_relax 0.895 [0.798, 0.967]) | — | Not in the ablation script output; to be added |
| Zenodo archive with DOI | — | To be added |

**Known housekeeping items.** Manifests record absolute paths from the authors' machine (`/Users/.../AnCore_v1.2/...`); Stage 6 reads receptor paths from the Stage 1 manifest, so re-running on another machine requires regenerating the workspace from Stage 1. The CLI silently ignores unrecognized flags (`parse_known_args`).

---

## 7. Installation and quick start

```bash
git clone https://github.com/norseenhosameldeen/AnCore-Scientific-Studio.git
cd AnCore-Scientific-Studio
conda env create -f environment.yml
conda activate ancore
```

**Run the full pipeline for one benchmark target** (writes to `workspaces/1STP/`):

```bash
python3 src/ancore_custom_runner.py --pdb-id 1STP --out-dir workspaces
```

**Run all five benchmark targets:**

```bash
python3 src/ancore_custom_runner.py --run-all --out-dir workspaces
```

> Flags accepted by the runner: `--pdb-id`, `--run-all`, `--exhaustiveness`, `--out-dir`, `--export-bundle`, `--bundle-only`, `--export-dir`. Other flags (for example `--stage`, `--project-dir`, `--pdb`, `--workspace`) are **not** recognized and are silently ignored; if no `--pdb-id` is given, the runner processes all five targets. Step-by-step execution and custom targets are provided through the desktop application.

**Re-compute ablations and controls from the shipped workspaces:**

```bash
python3 run_ablation_and_controls.py --workspaces-dir workspaces
python3 src/verify_q1_compliance.py --out-dir workspaces
```

**Build the macOS application** (requires [Rust](https://rustup.rs/) and `cargo install tauri-cli`):

```bash
cd ancore-tauri
./build_tauri_dmg.sh
```

---

## 8. Repository layout

```
├── src/                      # six-stage Python engine + integrity checker
├── ancore-tauri/             # Tauri/Rust shell and WKWebView UI (ui/index.html)
├── src-tauri/, ui/           # Tauri sources at repository root
├── workspaces/               # five-target benchmark outputs (manifests, CSV, PDB, figures, Mode C packages)
├── run_ablation_and_controls.py
├── environment.yml, requirements.txt
└── LICENSE                   # Apache-2.0
```

---

## 9. Citation

> Abdelhakim, N. H., & Senary, A. M. (2026). *AnCore: A Six-Stage Structure-Based Virtual Screening Pipeline Integrating Two-Phase Redocking, Reciprocal Rank Fusion Consensus Scoring, and a Contact-Anchored Elastic-Network Pose-Persistence Heuristic.* ChemRxiv. https://doi.org/10.26434/chemrxiv.15012010

```bibtex
@article{abdelhakim2026ancore,
  title   = {AnCore: A Six-Stage Structure-Based Virtual Screening Pipeline Integrating Two-Phase Redocking, Reciprocal Rank Fusion Consensus Scoring, and a Contact-Anchored Elastic-Network Pose-Persistence Heuristic},
  author  = {Abdelhakim, Norseen H. and Senary, Ahmed M.},
  year    = {2026},
  journal = {ChemRxiv},
  doi     = {10.26434/chemrxiv.15012010}
}
```

## 10. License

Distributed under the Apache License 2.0. See [`LICENSE`](LICENSE).
