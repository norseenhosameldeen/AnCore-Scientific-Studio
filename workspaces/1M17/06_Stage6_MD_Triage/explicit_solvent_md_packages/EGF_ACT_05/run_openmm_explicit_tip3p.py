#!/usr/bin/env python3
"""
AnCore v1.2 Mode C — Explicit-Solvent (TIP3P + 0.15 M NaCl) OpenMM Production MD Script
Target : 1M17
Ligand : Osimertinib (EGF_ACT_05) (EGF_ACT_05)
Force Fields: AMBER14 (amber14-all.xml + amber14/tip3pfb.xml) + OpenFF 2.0.0 (Sage) with AM1-BCC charges
"""
from pathlib import Path
import openmm as mm
from openmm import app, unit
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import SystemGenerator

def main():
    here = Path(__file__).resolve().parent
    receptor_pdb = Path("/Users/ahmedsenary/Desktop/AnCore_project/AnCore_v1.2/workspaces/1M17/01_Stage1_Target_and_Pocket/receptor_prepared_pH74_minimized.pdb")
    ligand_sdf = here / "EGF_ACT_05_docked.sdf"

    print("[Mode C] Loading ligand SDF and assigning AM1-BCC quantum-mechanical partial charges via OpenFF 2.0.0 (Sage)...")
    off_mol = Molecule.from_file(str(ligand_sdf))
    off_mol.assign_partial_charges(partial_charge_method="am1bcc")

    pdb = app.PDBFile(str(receptor_pdb))
    modeller = app.Modeller(pdb.topology, pdb.positions)
    lig_top = off_mol.to_topology().to_openmm()
    lig_pos = off_mol.conformers[0].to_openmm()
    modeller.add(lig_top, lig_pos)

    system_generator = SystemGenerator(
        forcefields=["amber14-all.xml", "amber14/tip3pfb.xml"],
        small_molecule_forcefield="openff-2.0.0",
        molecules=[off_mol],
        forcefield_kwargs={"constraints": app.HBonds, "rigidWater": True, "removeCMMotion": True},
        periodic_forcefield_kwargs={"nonbondedMethod": app.PME, "nonbondedCutoff": 1.0 * unit.nanometer},
    )

    print("[Mode C] Solvating in 10.0 A padded cubic TIP3P water box with 0.15 M NaCl physiological ionic strength...")
    modeller.addSolvent(
        system_generator.forcefield,
        model="tip3p",
        padding=1.0 * unit.nanometer,
        ionicStrength=0.15 * unit.molar,
    )

    system = system_generator.create_system(modeller.topology, molecules=[off_mol])
    system.addForce(mm.MonteCarloBarostat(1.0 * unit.bar, 300.0 * unit.kelvin, 25))

    integrator = mm.LangevinMiddleIntegrator(300.0 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds)
    sim = app.Simulation(modeller.topology, system, integrator)
    sim.context.setPositions(modeller.positions)

    print("[Mode C] Running 10,000-step L-BFGS energy minimization...")
    sim.minimizeEnergy(maxIterations=10000)

    print("[Mode C] Equilibrating (250 ps NPT at 300 K, 1.0 bar)...")
    sim.context.setVelocitiesToTemperature(300.0 * unit.kelvin)
    sim.step(125000)

    sim.reporters.append(app.DCDReporter(str(here / "1M17_EGF_ACT_05_production_10ns.dcd"), 5000))
    sim.reporters.append(app.StateDataReporter(
        str(here / "1M17_EGF_ACT_05_thermodynamics.csv"), 5000,
        step=True, time=True, potentialEnergy=True, temperature=True, density=True, speed=True
    ))
    print("[Mode C] Launching 10 ns explicit-solvent NPT production run (5,000,000 steps)...")
    sim.step(5000000)
    print("[Mode C] Production complete.")

if __name__ == "__main__":
    main()
