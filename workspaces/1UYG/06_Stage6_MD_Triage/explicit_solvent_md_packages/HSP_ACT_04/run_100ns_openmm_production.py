#!/usr/bin/env python3
"""
======================================================================================
AnCore v1.2 Mode C — Turnkey 100.0 ns (50,000,000 steps) OpenMM Production MD Launcher
Target : 1UYG
Ligand : Luminespib (HSP_ACT_04) (HSP_ACT_04)
Pre-Built Files Loaded Out-of-the-Box (Zero External Generator Dependencies Required):
  - System : openmm_100ns_system.xml (3332 full-atom particles)
  - State  : openmm_100ns_minimized_state.xml (E_min = -6958.09 kcal/mol)
  - PDB    : complex_100ns_full_atom.pdb
======================================================================================
"""
import argparse
from pathlib import Path
import openmm as mm
from openmm import app, unit

def main():
    parser = argparse.ArgumentParser(description="Run 100 ns OpenMM Production MD for 1UYG + HSP_ACT_04")
    parser.add_argument("--ns", type=float, default=100.0, help="Production duration in nanoseconds (default: 100.0 ns)")
    parser.add_argument("--report-ps", type=float, default=20.0, help="Reporter interval in picoseconds (default: 20.0 ps)")
    args = parser.parse_args()

    here = Path(__file__).resolve().parent
    pdb_file = app.PDBFile(str(here / "complex_100ns_full_atom.pdb"))
    system = mm.XmlSerializer.deserialize((here / "openmm_100ns_system.xml").read_text(encoding="utf-8"))
    state = mm.XmlSerializer.deserialize((here / "openmm_100ns_minimized_state.xml").read_text(encoding="utf-8"))

    dt_ps = 0.002  # 2.0 fs timestep
    total_steps = int(round((args.ns * 1000.0) / dt_ps))
    report_steps = max(100, int(round(args.report_ps / dt_ps)))

    integrator = mm.LangevinMiddleIntegrator(300.0 * unit.kelvin, 1.0 / unit.picosecond, dt_ps * unit.picoseconds)
    sim = app.Simulation(pdb_file.topology, system, integrator)
    sim.context.setState(state)

    dcd_out = here / f"1UYG_HSP_ACT_04_production_{int(args.ns)}ns.dcd"
    csv_out = here / f"1UYG_HSP_ACT_04_thermodynamics_{int(args.ns)}ns.csv"
    chk_out = here / f"1UYG_HSP_ACT_04_checkpoint.chk"

    sim.reporters.append(app.DCDReporter(str(dcd_out), report_steps))
    sim.reporters.append(app.StateDataReporter(
        str(csv_out), report_steps,
        step=True, time=True, potentialEnergy=True, kineticEnergy=True,
        totalEnergy=True, temperature=True, speed=True, progress=True,
        remainingTime=True, totalSteps=total_steps
    ))
    sim.reporters.append(app.CheckpointReporter(str(chk_out), report_steps * 10))

    print(f"[AnCore Mode C] Launching {args.ns:.1f} ns ({total_steps:,} steps at 2.0 fs) production MD for 1UYG + HSP_ACT_04...")
    sim.step(total_steps)
    final_pdb = here / f"1UYG_HSP_ACT_04_final_{int(args.ns)}ns.pdb"
    st_end = sim.context.getState(getPositions=True)
    with open(final_pdb, "w", encoding="utf-8") as fh:
        app.PDBFile.writeFile(pdb_file.topology, st_end.getPositions(), fh)
    print(f"[AnCore Mode C] Completed {args.ns:.1f} ns run -> {dcd_out}")

if __name__ == "__main__":
    main()
