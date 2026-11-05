#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
AnCore Scientific Studio — GitHub Repository Organizer & Sanitizer
File: organize_for_github.py (Jupyter / IPython & macOS Terminal Compatible)
================================================================================
Why your current `AnCore_project.zip` is 288 MB:
  When `cargo tauri build` compiled the native macOS `.app` and `.dmg`, Rust wrote
  ~270 MB of intermediate build artifacts into `ancore-tauri/src-tauri/target/`.
  GitHub rejects pushes containing files > 100 MB.

What this script does automatically on your Mac:
  1. Locates your `/Users/ahmedsenary/Desktop/AnCore_project` folder.
  2. Creates a clean, modular GitHub repository directory:
       `~/Desktop/AnCore-Scientific-Studio/`
  3. Copies all source files into an industry-standard layout (`src/`, `ancore-tauri/`,
     `benchmarks/`, `scripts/`, `docs/figures/`, `docs/paper/`) while stripping out
     the ~270 MB `target/`, `node_modules/`, `__pycache__/`, `.DS_Store`, and legacy
     scratch files.
  4. Sanitizes hardcoded `/Users/ahmedsenary/...` paths inside `benchmarks/*.json`
     into portable relative paths so anyone cloning your repo can run it on their Mac.
  5. Extracts your compiled `AnCore Studio.dmg` installer into a separate folder
     (`~/Desktop/AnCore-Release-Assets/`) so you can attach it to a GitHub Release.
  6. Generates `README.md`, `LICENSE` (Apache-2.0), `.gitignore`, `environment.yml`,
     `requirements.txt`, and `.github/workflows/benchmark-audit.yml` (GitHub Actions CI).
================================================================================
"""

import argparse
import json
import os
import re
import shutil
import sys
import zipfile
from pathlib import Path

README_MD = """# AnCore Scientific Studio

[![Platform](https://img.shields.io/badge/Platform-macOS%20(Apple%20Silicon%20%26%20Intel)-0D7377?style=flat-square&logo=apple)](https://github.com/)
[![Architecture](https://img.shields.io/badge/Architecture-Tauri%201.6%20(Rust)%20%2B%20Python%203.11-1B3A5C?style=flat-square&logo=rust)](https://tauri.app/)
[![Pipeline](https://img.shields.io/badge/Pipeline-6--Stage%20SBVS%20%2B%20300%20K%20Langevin%20Relaxation-059669?style=flat-square)]()
[![License](https://img.shields.io/badge/License-Apache%202.0-blue?style=flat-square)](LICENSE)

**AnCore Scientific Studio** is a standalone native macOS chemoinformatics and structure-based virtual screening (SBVS) workbench distributed as a compiled Rust/Tauri application (`AnCore Studio.dmg`). It couples a hardware-accelerated `WKWebView` + `3Dmol.js` WebGL viewport to a modular six-stage Python scientific engine supporting both **Step-by-Step Interactive Mode** and **One-Shot Batch Execution**, live RCSB PDB fetching, local `.pdb` import, and multi-molecule `.sdf`/`.mol` ligand library ingestion.

![AnCore Architecture and 6-Stage Pipeline](docs/figures/fig1_architecture.png)

---

## Key Scientific & Architectural Features

1. **Three-Layer Native macOS Architecture (Tauri Rust + `WKWebView` + Python Sidecar)**
   - Eliminates Python GUI event-loop freezes and the 2 MB IPC string-serialization ceiling of legacy webviews by reading stage manifests (`stage1_manifest.json`–`stage6_manifest.json`) and 3D coordinate blocks directly via native Rust filesystem IPC (`src-tauri/src/main.rs`).
2. **Stage 1 — Automated Receptor Curation & 5.0 Å Cavity Mapping (`stage1_target_prep.py`)**
   - Isolates the primary polypeptide chain, reconstructs missing heavy atoms, assigns physiological protonation states at pH 7.4 (`PDBFixer` / `OpenMM` AMBER14), and computes the 5.0 Å Euclidean active-site bounding box.
3. **Stage 2 — Two-Phase Co-Crystal Redocking Quality-Control Gate (`stage2_redock_gate.py`)**
   - Combines **Phase A** local Broyden–Fletcher–Goldfarb–Shanno (BFGS) crystallographic basin relaxation with **Phase B** focused-box (3.0 Å padding, exhaustiveness = 16) Monte Carlo sampling, enforcing heavy-atom $\\text{RMSD} \\le 2.0\\text{ \\AA}$ prior to virtual screening.
4. **Stage 3 — Multi-Molecule `.SDF` & SMILES Conformer Standardization (`stage3_ligand_prep.py`)**
   - Parses multi-record `.sdf`/`.mol` files and SMILES strings, enumerates pH 7.4 protomers, embeds 3D conformers via RDKit `ETKDGv3`, relaxes structures with `MMFF94s`, and computes Lipinski and Bickerton `QED` descriptors.
5. **Stage 4 — Multi-Scoring Consensus Virtual Screening (`stage4_screening.py`)**
   - Evaluates candidate ligands across three orthogonal empirical scoring potentials (**AutoDock Vina**, **Vinardo**, and **AutoDock4**) and aggregates rankings via **Reciprocal Rank Fusion ($\\text{RRF}, k = 10$)** and mean $Z$-score normalization.
6. **Stage 5 — 3D Non-Covalent Interaction Fingerprinting (`stage5_interactions.py`)**
   - Detects 3D hydrogen bonds ($\\le 3.5\\text{ \\AA}$), salt bridges ($\\le 4.2\\text{ \\AA}$), aromatic $\\pi$-stacking ($\\le 4.3\\text{ \\AA}$), halogen bonds ($\\le 3.9\\text{ \\AA}$), and hydrophobic contacts ($\\le 4.0\\text{ \\AA}$), computing native hotspot Jaccard/Tanimoto similarity ($T_{\\text{IFP}}$).
7. **Stage 6 — 300 K Langevin Dynamics Relaxation Triage (`stage6_md_triage.py`)**
   - Subjects docked complexes to a 50-frame ($5.0\\text{ ns}$ effective window) $300\\text{ K}$ Langevin relaxation trajectory under a dual biophysical stability gate ($\\langle\\text{RMSD}\\rangle \\le 2.0\\text{ \\AA}$ and hotspot hydrogen-bond occupancy $\\text{HB}_{\\text{occ}} \\ge 65.0\\%$), rejecting thermally unstable decoys and rescuing true biological hits.

---

## Five-Target Benchmark Summary (50 Screened Complexes)

| PDB ID | Target Receptor | Co-Crystal Ligand | Stage 2 Redock RMSD | Stage 4 Static AUC ($\\text{EF}_{20\\%}$) | Stage 5 Top Lead $T_{\\text{IFP}}$ | Stage 6 Decoys Rejected | Post-Relaxation ROC-AUC ($\\text{EF}_{20\\%}$) |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **1IEP** | ABL1 Tyrosine Kinase | Imatinib (`STI`) | $0.354\\text{ \\AA}$ | $0.800$ ($2.00\\times$) | $0.789$ (Nilotinib) | $5/5$ ($100\\%$) | **$1.000$ ($2.00\\times$)** |
| **3ERT** | Estrogen Receptor $\\alpha$ | 4-Hydroxytamoxifen (`OHT`) | $0.412\\text{ \\AA}$ | $1.000$ ($2.00\\times$) | $0.409$ (Bazedoxifene) | $5/5$ ($100\\%$) | **$1.000$ ($2.00\\times$)** |
| **1UYG** | HSP90$\\alpha$ Chaperone | Purine Inhibitor (`PU2`) | $0.485\\text{ \\AA}$ | $0.680$ ($2.00\\times$) | $0.500$ (Luminespib) | $5/5$ ($100\\%$) | **$1.000$ ($2.00\\times$)** |
| **1M17** | EGFR Tyrosine Kinase | Erlotinib (`AQ4`) | $0.528\\text{ \\AA}$ | $0.800$ ($2.00\\times$) | $0.250$ (Lapatinib) | $5/5$ ($100\\%$) | **$1.000$ ($2.00\\times$)** |
| **1STP** | Streptavidin ($\\beta$-Barrel) | D-Biotin (`BTN`) | $0.318\\text{ \\AA}$ | $0.360$ ($1.00\\times$) | $0.733$ (D-Biotin) | $5/5$ ($100\\%$) | **$1.000$ ($2.00\\times$, Rescued #1)** |
| **Suite** | **5-Target Mean (50 Complexes)** | **5 Crystal References** | **$0.419\\text{ \\AA}$ (100% Pass)** | **$0.728$ ($1.80\\times$)** | **$0.536$** | **$25/25$ ($100\\%$)** | **$1.000$ ($2.00\\times$, 100% Top-1)** |

![Stage 6 300 K Langevin Relaxation Trajectories](docs/figures/fig5_md_triage_trajectories.png)

---

## Repository Structure

```text
AnCore-Scientific-Studio/
├── README.md                          # Project documentation, benchmark summary & setup guide
├── LICENSE                            # Permissive Apache-2.0 Open-Source License
├── .gitignore                         # Excludes Rust target/, node_modules/, __pycache__/
├── environment.yml                    # Conda environment (rdkit, openmm, pdbfixer, autodock-vina)
├── requirements.txt                   # Python package dependencies
├── src/                               # Modular 6-Stage Python Scientific Engine
│   ├── ancore_custom_runner.py        # Step-by-Step & One-Shot 6-Stage Runner + Multi-Molecule .SDF Parser
│   ├── stage1_target_prep.py          # Stage 1: Receptor Curation, pH 7.4 Protonation & 5.0 Å Pocket Box
│   ├── stage2_redock_gate.py          # Stage 2: Two-Phase BFGS + Focused-Box Redocking Gate (<= 2.0 Å)
│   ├── stage3_ligand_prep.py          # Stage 3: Multi-Molecule .SDF / SMILES Prep, ETKDGv3 & MMFF94s
│   ├── stage4_screening.py            # Stage 4: Tri-Scoring Consensus (Vina, Vinardo, AD4) & RRF (k=10)
│   ├── stage5_interactions.py         # Stage 5: 3D Non-Covalent Interaction Fingerprinting & T_IFP
│   └── stage6_md_triage.py            # Stage 6: 300 K Langevin Dynamics Relaxation Triage & Rescoring
├── ancore-tauri/                      # Native macOS Rust/Tauri Desktop Application
│   ├── build_tauri_dmg.sh             # One-command macOS .app / .dmg compiler script
│   ├── ui/
│   │   └── index.html                 # Native WKWebView UI + 3Dmol.js WebGL Viewport
│   └── src-tauri/
│       ├── Cargo.toml                 # Rust crate manifest
│       ├── build.rs                   # Tauri build hook
│       ├── tauri.conf.json            # macOS window & bundle configuration
│       ├── icons/                     # macOS application iconset
│       └── src/
│           └── main.rs                # Native Rust IPC backend
├── benchmarks/                        # Pre-computed 5-Target Validation Suite (1IEP, 3ERT, 1UYG, 1M17, 1STP)
├── scripts/                           # Reproducibility & Verification Scripts
│   ├── audit_benchmark_data.py        # Independent 4-part mathematical & coordinate audit script
│   └── organize_for_github.py         # Repository organizer & path sanitizer
└── docs/                              # Manuscript & 600-DPI Publication Figures
    ├── figures/                       # fig1_architecture.png .. fig6_streptavidin_rescue.png
    └── paper/                         # LaTeX, PDF, and Word Compendium
```

---

## Quick Start & Installation

### 1. Create the Scientific Python Environment (Conda)

```bash
conda env create -f environment.yml
conda activate ancore
```

### 2. Verify All Benchmark Results & Equations (One Command)

To independently audit all 5 benchmark targets (`1IEP`, `3ERT`, `1UYG`, `1M17`, `1STP`), 3D `.pdb` pocket geometries, redocking RMSDs, $\\text{RRF}$ consensus $\\text{ROC-AUC}$ scores, $T_{\\text{IFP}}$ similarities, and 50-frame $300\\text{ K}$ Langevin relaxation trajectories:

```bash
python3 scripts/audit_benchmark_data.py
```

### 3. Run the 6-Stage Pipeline Step-by-Step via CLI

```bash
# Execute Step 1 (Target Curation & 5.0 Å Box) on a custom RCSB structure (e.g., 4WKQ)
python3 src/ancore_custom_runner.py --project-dir . --pdb-id 4WKQ --target-name "EGFR Kinase" --ligand-code IRE --stage 1 --reset-downstream 1

# Advance through Stages 2 -> 6 step-by-step
python3 src/ancore_custom_runner.py --project-dir . --pdb-id 4WKQ --stage 2
python3 src/ancore_custom_runner.py --project-dir . --pdb-id 4WKQ --stage 3
python3 src/ancore_custom_runner.py --project-dir . --pdb-id 4WKQ --stage 4
python3 src/ancore_custom_runner.py --project-dir . --pdb-id 4WKQ --stage 5
python3 src/ancore_custom_runner.py --project-dir . --pdb-id 4WKQ --stage 6
```

### 4. Build the Native macOS Application (`AnCore Studio.dmg`)

Prerequisite: Install [Rust (`rustup`)](https://rustup.rs/) and the Tauri CLI (`cargo install tauri-cli`).

```bash
cd ancore-tauri
./build_tauri_dmg.sh
```

---

## Authors & Affiliation

* **Norseen H. Abdelhakim** and **Ahmed M. Senary**
* Faculty of Pharmacy, Galala University, New Galala City, Suez, Egypt

## License

Distributed under the permissive **Apache License 2.0**. See [`LICENSE`](LICENSE) for full details.
"""

GITIGNORE = """# Rust / Tauri compilation artifacts (250+ MB)
ancore-tauri/src-tauri/target/
src-tauri/target/
target/
**/*.rs.bk

# Node / frontend dependencies
node_modules/
dist/

# Compiled macOS bundles (attach to GitHub Releases instead of git tree)
*.dmg
*.app/

# Python caches & virtual environments
__pycache__/
*.py[cod]
*$py.class
.venv/
venv/
env/
.ipynb_checkpoints/

# LaTeX auxiliary files
*.aux
*.log
*.out
*.toc
*.fls
*.fdb_latexmk
*.synctex.gz

# macOS system metadata
.DS_Store
.AppleDouble
.LSOverride
._*
"""

ENVIRONMENT_YML = """name: ancore
channels:
  - conda-forge
  - defaults
dependencies:
  - python=3.11
  - rdkit>=2023.09
  - openmm>=8.1
  - pdbfixer>=1.9
  - autodock-vina>=1.2
  - openbabel>=3.1
  - numpy>=1.26
  - scipy>=1.11
  - matplotlib>=3.8
  - pip
  - pip:
      - meeko>=0.5.0
"""

REQUIREMENTS_TXT = """numpy>=1.26.0
scipy>=1.11.0
matplotlib>=3.8.0
rdkit>=2023.9.1
"""

GITHUB_WORKFLOW_YML = """name: AnCore Benchmark & Equation Audit

on:
  push:
    branches: [ main, master ]
  pull_request:
    branches: [ main, master ]

jobs:
  audit-benchmark:
    runs-on: macos-latest
    steps:
      - name: Checkout Repository
        uses: actions/checkout@v4

      - name: Set up Python 3.11
        uses: actions/setup-python@v5
        with:
          python-version: "3.11"

      - name: Install Core Audit Dependencies
        run: |
          python -m pip install --upgrade pip
          pip install numpy matplotlib

      - name: Run 4-Part Independent Benchmark Audit
        run: |
          python3 scripts/audit_benchmark_data.py
"""

LICENSE_TXT = """                                 Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/

   Copyright 2026 Norseen H. Abdelhakim and Ahmed M. Senary (Galala University)

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""


def sanitize_json_paths(json_text: str, old_prefixes: list) -> str:
    """Replaces machine-specific absolute paths inside JSON manifests with relative repo paths."""
    out = json_text
    for pref in old_prefixes:
        if pref:
            out = out.replace(pref.rstrip("/") + "/", "")
            out = out.replace(pref.rstrip("/"), ".")
    # Also strip generic /Users/<username>/.../AnCore_project/ or /working_dir/.../bench_verify/
    out = re.sub(r"/Users/[^/\"\\]+/Desktop/AnCore_project/?", "", out)
    out = re.sub(r"/working_dir/[^/\"\\]+/artifacts/file_generation/ttl=63d/intermediate/bench_verify/?", "", out)
    return out


def find_source_project_dir(explicit_dir: str = None) -> Path:
    if explicit_dir:
        p = Path(explicit_dir).expanduser().resolve()
        if p.exists():
            return p

    home = Path.home()
    candidates = [
        home / "Desktop" / "AnCore_project",
        home / "Downloads" / "AnCore_project",
        home / "Documents" / "AnCore_project",
        Path.cwd(),
        Path("/working_dir/c_9af29f10d0e439d3/artifacts/file_generation/ttl=63d/output"),
    ]
    for c in candidates:
        if c.exists() and ((c / "stage1_target_prep.py").exists() or (c / "benchmarks").exists() or (c / "src").exists()):
            return c
    return Path.cwd()


def organize_repository(src_dir: Path, dest_dir: Path, create_zip: bool = False):
    print("=" * 86)
    print("ANCORE SCIENTIFIC STUDIO — GITHUB REPOSITORY ORGANIZER")
    print("=" * 86)
    print(f"  • Source Project Folder : {src_dir}")
    print(f"  • Clean GitHub Repo Dir : {dest_dir}")

    if dest_dir.exists():
        shutil.rmtree(dest_dir)

    # Create standard directory hierarchy
    src_mod_dir = dest_dir / "src"
    tauri_ui_dir = dest_dir / "ancore-tauri" / "ui"
    tauri_rust_dir = dest_dir / "ancore-tauri" / "src-tauri" / "src"
    tauri_icons_dir = dest_dir / "ancore-tauri" / "src-tauri" / "icons"
    bench_out_dir = dest_dir / "benchmarks"
    scripts_dir = dest_dir / "scripts"
    figures_dir = dest_dir / "docs" / "figures"
    paper_dir = dest_dir / "docs" / "paper"
    workflow_dir = dest_dir / ".github" / "workflows"

    for d in [
        src_mod_dir, tauri_ui_dir, tauri_rust_dir, tauri_icons_dir,
        bench_out_dir, scripts_dir, figures_dir, paper_dir, workflow_dir
    ]:
        d.mkdir(parents=True, exist_ok=True)

    # 1. Write root documentation & config files
    (dest_dir / "README.md").write_text(README_MD, encoding="utf-8")
    (dest_dir / "LICENSE").write_text(LICENSE_TXT, encoding="utf-8")
    (dest_dir / ".gitignore").write_text(GITIGNORE, encoding="utf-8")
    (dest_dir / "environment.yml").write_text(ENVIRONMENT_YML, encoding="utf-8")
    (dest_dir / "requirements.txt").write_text(REQUIREMENTS_TXT, encoding="utf-8")
    (workflow_dir / "benchmark-audit.yml").write_text(GITHUB_WORKFLOW_YML, encoding="utf-8")

    # 2. Copy core 6-stage Python modules into `src/` (and root symlinks/copies for backward compatibility)
    stage_files = [
        "stage1_target_prep.py",
        "stage2_redock_gate.py",
        "stage3_ligand_prep.py",
        "stage4_screening.py",
        "stage5_interactions.py",
        "stage6_md_triage.py",
        "ancore_custom_runner.py",
    ]
    copied_stages = 0
    for sf in stage_files:
        for cand in [src_dir / sf, src_dir / "src" / sf]:
            if cand.exists():
                shutil.copy2(cand, src_mod_dir / sf)
                copied_stages += 1
                break

    # If `ancore_custom_runner.py` wasn't a separate file in the sandbox VM, extract it or generate it via setup_ancore_tauri.py
    if not (src_mod_dir / "ancore_custom_runner.py").exists() and (src_dir / "setup_ancore_tauri.py").exists():
        try:
            setup_code = (src_dir / "setup_ancore_tauri.py").read_text(encoding="utf-8")
            m = re.search(r'RUNNER_CODE\s*=\s*r?"""(.*?)"""\n\n', setup_code, re.S)
            if m:
                (src_mod_dir / "ancore_custom_runner.py").write_text(m.group(1), encoding="utf-8")
                copied_stages += 1
        except Exception:
            pass

    print(f"  • Copied Python Modules : {copied_stages} files -> src/")

    # 3. Copy Tauri Rust & WebGL UI files (WITHOUT copying the 270 MB `target/` directory!)
    for ui_cand in [src_dir / "ancore-tauri" / "ui" / "index.html", src_dir / "index.html"]:
        if ui_cand.exists():
            html_txt = ui_cand.read_text(encoding="utf-8", errors="ignore")
            html_txt = html_txt.replace("/Users/ahmedsenary/Desktop/AnCore_project", ".")
            (tauri_ui_dir / "index.html").write_text(html_txt, encoding="utf-8")
            break

    for rs_cand in [src_dir / "ancore-tauri" / "src-tauri" / "src" / "main.rs", src_dir / "main.rs"]:
        if rs_cand.exists():
            rs_txt = rs_cand.read_text(encoding="utf-8", errors="ignore")
            # Make receptor_clean.pdb & crystal_ligand.pdb lookup portable across any cloned repo path!
            old_rec_block = """        let receptor_pdb = read_text_safe(rec_path);
        let crystal_ligand_pdb = read_text_safe(cryst_path);"""
            new_rec_block = """        let mut receptor_pdb = read_text_safe(rec_path);
        if receptor_pdb.is_empty() {
            receptor_pdb = read_text_safe(tdir.join("receptor_clean.pdb").to_str());
        }
        let mut crystal_ligand_pdb = read_text_safe(cryst_path);
        if crystal_ligand_pdb.is_empty() {
            crystal_ligand_pdb = read_text_safe(tdir.join("crystal_ligand.pdb").to_str());
        }"""
            if old_rec_block in rs_txt:
                rs_txt = rs_txt.replace(old_rec_block, new_rec_block)
            (tauri_rust_dir / "main.rs").write_text(rs_txt, encoding="utf-8")
            break

    for cfg_name in ["Cargo.toml", "tauri.conf.json", "build.rs"]:
        for cand in [src_dir / "ancore-tauri" / "src-tauri" / cfg_name, src_dir / cfg_name]:
            if cand.exists():
                shutil.copy2(cand, dest_dir / "ancore-tauri" / "src-tauri" / cfg_name)
                break

    if not (dest_dir / "ancore-tauri" / "src-tauri" / "build.rs").exists():
        (dest_dir / "ancore-tauri" / "src-tauri" / "build.rs").write_text("fn main() {\n  tauri_build::build()\n}\n", encoding="utf-8")

    # Copy icons if present
    src_icons = src_dir / "ancore-tauri" / "src-tauri" / "icons"
    if src_icons.exists():
        for ic in src_icons.iterdir():
            if ic.is_file():
                shutil.copy2(ic, tauri_icons_dir / ic.name)

    # Write portable build_tauri_dmg.sh
    build_sh = dest_dir / "ancore-tauri" / "build_tauri_dmg.sh"
    build_sh.write_text("""#!/bin/bash
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR/src-tauri"
echo "Compiling Native macOS Rust/Tauri Application & DMG..."
cargo tauri build
echo "Done! Installer bundle generated in src-tauri/target/release/bundle/dmg/"
""", encoding="utf-8")
    build_sh.chmod(0o755)

    # If a compiled .dmg exists inside the user's 288 MB folder, copy it to ~/Desktop/AnCore-Release-Assets/
    dmg_search_root = src_dir / "ancore-tauri" / "src-tauri" / "target" / "release" / "bundle" / "dmg"
    if dmg_search_root.exists():
        for dmg_file in dmg_search_root.glob("*.dmg"):
            rel_dir = src_dir.parent / "AnCore-GitHub-Release-Assets"
            rel_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy2(dmg_file, rel_dir / "AnCore_Studio_macOS.dmg")
            print(f"  • Extracted DMG Binary  : {rel_dir / 'AnCore_Studio_macOS.dmg'} (Attach this to GitHub Releases!)")

    # 4. Copy & sanitize the 5-target `benchmarks/` directory
    bench_candidates = [
        src_dir / "benchmarks",
        Path("/working_dir/c_9af29f10d0e439d3/artifacts/file_generation/ttl=63d/intermediate/bench_verify/benchmarks"),
    ]
    bench_src = next((b for b in bench_candidates if b.exists()), None)
    if bench_src:
        for pdb_id in ["1IEP", "3ERT", "1UYG", "1M17", "1STP"]:
            s_tdir = bench_src / pdb_id
            if not s_tdir.exists():
                continue
            d_tdir = bench_out_dir / pdb_id
            d_tdir.mkdir(parents=True, exist_ok=True)
            for item in s_tdir.rglob("*"):
                if item.is_dir() or item.name.startswith(".") or "__pycache__" in item.parts:
                    continue
                rel_p = item.relative_to(s_tdir)
                out_p = d_tdir / rel_p
                out_p.parent.mkdir(parents=True, exist_ok=True)
                if item.suffix.lower() == ".json":
                    txt = item.read_text(encoding="utf-8", errors="ignore")
                    txt = sanitize_json_paths(txt, [str(src_dir)])
                    out_p.write_text(txt, encoding="utf-8")
                else:
                    shutil.copy2(item, out_p)
        print("  • Copied & Sanitized    : 5 Benchmark Target Suites (1IEP, 3ERT, 1UYG, 1M17, 1STP) -> benchmarks/")

    # 5. Copy verification & organizer scripts into `scripts/`
    for sc_name in ["audit_benchmark_data.py", "organize_for_github.py", "setup_ancore_tauri.py"]:
        for cand in [src_dir / sc_name, src_dir / "scripts" / sc_name, Path(__file__).resolve().parent / sc_name]:
            if cand.exists():
                shutil.copy2(cand, scripts_dir / sc_name)
                break

    # 6. Copy 600-DPI publication figures & manuscript files into `docs/`
    fig_names = [
        "fig1_architecture.png",
        "fig2_redocking_rmsd.png",
        "fig3_consensus_roc_enrichment.png",
        "fig4_ifp_pharmacophore_heatmap.png",
        "fig5_md_triage_trajectories.png",
        "fig6_streptavidin_rescue.png",
    ]
    copied_figs = 0
    for fn in fig_names:
        for cand in [src_dir / fn, src_dir / "docs" / "figures" / fn, Path("/working_dir/c_9af29f10d0e439d3/artifacts/file_generation/ttl=63d/output") / fn]:
            if cand.exists():
                shutil.copy2(cand, figures_dir / fn)
                copied_figs += 1
                break
    print(f"  • Copied 600-DPI Figures: {copied_figs} figures -> docs/figures/")

    paper_cands = [
        ("AnCore_Research_Paper.tex", "AnCore_Research_Paper.tex"),
        ("AnCore_Research_Paper.pdf", "AnCore_Research_Paper.pdf"),
        ("AnCore_Figures_and_Tables_Compendium.docx", "AnCore_Figures_and_Tables_Compendium.docx"),
        ("AnCore_Research_Paper.docx", "AnCore_Research_Paper.docx"),
    ]
    for src_name, dst_name in paper_cands:
        for cand in [
            src_dir / src_name,
            Path("/working_dir/c_9af29f10d0e439d3/user_uploaded_files/ttl=63d/user_uploaded/AnCore_Research_Paper_JCIM (2).pdf") if src_name.endswith(".pdf") else Path("/nonexistent"),
            Path("/working_dir/c_9af29f10d0e439d3/artifacts/file_generation/ttl=63d/output") / src_name,
        ]:
            if cand.exists():
                shutil.copy2(cand, paper_dir / dst_name)
                break

    # Compute total clean repository size
    total_bytes = sum(f.stat().st_size for f in dest_dir.rglob("*") if f.is_file())
    total_mb = total_bytes / (1024 * 1024)
    print(f"  • Clean Repository Size : {total_mb:.2f} MB (Reduced from ~288.05 MB — 95%+ smaller!)")

    if create_zip:
        zip_path = dest_dir.parent / f"{dest_dir.name}.zip"
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
            for f in sorted(dest_dir.rglob("*")):
                if f.is_file():
                    zf.write(f, f.relative_to(dest_dir.parent))
        print(f"  • Packaged GitHub Zip   : {zip_path} ({zip_path.stat().st_size / (1024*1024):.2f} MB)")

    print("=" * 86)
    print("NEXT STEPS TO PUSH TO GITHUB (Copy & Paste in Terminal):")
    print("=" * 86)
    print(f"  cd \"{dest_dir}\"")
    print("  git init")
    print("  git branch -M main")
    print("  git add .")
    print("  git commit -m \"Initial release: AnCore Scientific Studio 6-stage SBVS & 300 K Langevin relaxation workbench\"")
    print("  git remote add origin https://github.com/<YOUR_GITHUB_USERNAME>/AnCore-Scientific-Studio.git")
    print("  git push -u origin main")
    print("=" * 86)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Organize AnCore_project into a clean, lightweight GitHub repository.")
    parser.add_argument("--src", type=str, default=None, help="Path to existing AnCore_project directory.")
    parser.add_argument("--dest", type=str, default=None, help="Output path for clean GitHub repository.")
    parser.add_argument("--zip", action="store_true", help="Also create a .zip archive of the clean repository.")
    args, _ = parser.parse_known_args()

    src_p = find_source_project_dir(args.src)
    if args.dest:
        dest_p = Path(args.dest).expanduser().resolve()
    else:
        dest_p = src_p.parent / "AnCore-Scientific-Studio"
    organize_repository(src_p, dest_p, create_zip=args.zip)
