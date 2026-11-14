#!/usr/bin/env bash
set -euo pipefail
# AnCore v1.2 Mode C — Turnkey 100.0 ns (50,000,000 steps) GROMACS 2024 Production Pipeline
# Target: 3ERT | Top Hit: Raloxifene (ERA_ACT_02) (ERA_ACT_02)
echo "[AnCore Mode C] Preparing explicit TIP3P solvated box and 0.15 M NaCl ions for 3ERT + ERA_ACT_02..."
gmx solvate -cp complex_100ns_ready.gro -cs spc216.gro -p topol.top -o solvated.gro
gmx grompp -f gromacs_100ns_npt.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 2
echo "SOL" | gmx genion -s ions.tpr -o solvated_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15
gmx grompp -f gromacs_100ns_npt.mdp -c solvated_ions.gro -p topol.top -o prod_100ns.tpr -maxwarn 2
echo "[AnCore Mode C] Launching 100.0 ns NPT production MD (50,000,000 steps)..."
gmx mdrun -deffnm prod_100ns -v
