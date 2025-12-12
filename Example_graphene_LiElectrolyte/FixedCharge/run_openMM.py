from __future__ import print_function
import os
import sys
from os import path

sys.path.append('../../lib/')
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from MM_classes_FV import *
from Fixed_Voltage_routines import *

# basic run control
simulation_time_ns = float(os.environ.get("SIM_TIME_NS", 1.0))
freq_charge_update_fs = int(os.environ.get("FREQ_CHARGE_UPDATE_FS", 50))
freq_traj_output_ps = int(os.environ.get("FREQ_TRAJ_OUTPUT_PS", 50))
freq_checkpoint_ps = int(os.environ.get("FREQ_CHECKPOINT_PS", 10))
write_charges = os.environ.get("WRITE_CHARGES", "true").lower() == "true"

checkpoint = 'state.chk'
charge_name = 'charges.dat'

# simulation mode and voltage
simulation_type = os.environ.get("SIM_MODE", "MC_equil")  # "MC_equil" or "Constant_V"
Voltage = float(os.environ.get("APPLIED_VOLTAGE", 0.0))

# electrode identifiers (chain indices)
cathode_index = (0, 2)
anode_index = (1, 3)

# input files
graphene_ffdir = os.environ.get(
    "GRAPHENE_FF_DIR",
    "../../Example_graphene_BMIM_BF4_ACN_10pct/graphene_ffdir/"
)
electrolyte_ff = os.environ.get("ELECTROLYTE_FF_XML", "../forcefield/li_fixed_charge.xml")
electrolyte_residues = os.environ.get("ELECTROLYTE_RESIDUE_XML", "../forcefield/li_fixed_charge_residues.xml")

packmol_pdb = os.environ.get("PACKMOL_PDB", "../packmol/li_electrolyte_box.pdb")
equilibrated_pdb = os.environ.get("EQUIL_PDB", "equilibrated.pdb")
input_pdb = equilibrated_pdb if path.exists(equilibrated_pdb) else packmol_pdb
input_pdb = os.environ.get("INPUT_PDB", input_pdb)

# build system
MMsys = MM_FixedVoltage(
    pdb_list=[input_pdb],
    residue_xml_list=[
        electrolyte_residues,
        graphene_ffdir + 'graph_residue_c.xml',
        graphene_ffdir + 'graph_residue_n.xml'
    ],
    ff_xml_list=[
        electrolyte_ff,
        graphene_ffdir + 'graph.xml',
        graphene_ffdir + 'graph_c_freeze.xml',
        graphene_ffdir + 'graph_n_freeze.xml'
    ]
)

MMsys.set_periodic_residue(True)
platform_name = os.environ.get("FV_PLATFORM", "Reference")
MMsys.set_platform(platform_name)

MMsys.initialize_electrodes(
    Voltage,
    cathode_identifier=cathode_index,
    anode_identifier=anode_index,
    chain=True,
    exclude_element=("H",)
)
MMsys.initialize_electrolyte(Natom_cutoff=100)
MMsys.generate_exclusions(flag_SAPT_FF_exclusions=False)

# allow restart
if path.exists(checkpoint):
    print("restarting simulation from checkpoint ", checkpoint)
    MMsys.simmd.loadCheckpoint(checkpoint)

state = MMsys.simmd.context.getState(getEnergy=True, getForces=True, getVelocities=False, getPositions=True)
positions = state.getPositions()
with open('start.pdb', 'w') as handle:
    PDBFile.writeFile(MMsys.simmd.topology, positions, handle)

append_trajectory = False
trajectory_file_name = 'equil_MC.dcd' if simulation_type == "MC_equil" else 'FV_NVT.dcd'
if simulation_type == "Constant_V" and path.exists(trajectory_file_name):
    append_trajectory = True

MMsys.set_trajectory_output(
    trajectory_file_name,
    freq_traj_output_ps * 1000,
    append_trajectory,
    checkpoint,
    freq_checkpoint_ps * 1000
)

if simulation_type == "MC_equil":
    celldim = MMsys.simmd.topology.getUnitCellDimensions()
    MMsys.MC = MC_parameters(
        MMsys.temperature,
        celldim,
        electrode_move="Anode",
        pressure=1.0 * bar,
        barofreq=100,
        shiftscale=0.2
    )

for i in range(int(simulation_time_ns * 1000 / freq_traj_output_ps)):
    state = MMsys.simmd.context.getState(getEnergy=True, getForces=True, getVelocities=False, getPositions=True)
    print(str(state.getKineticEnergy()))
    print(str(state.getPotentialEnergy()))
    for j in range(MMsys.system.getNumForces()):
        f = MMsys.system.getForce(j)
        print(type(f), str(MMsys.simmd.context.getState(getEnergy=True, groups=2 ** j).getPotentialEnergy()))
    print("", flush=True)

    if simulation_type == "MC_equil":
        for j in range(int(freq_traj_output_ps * 1000 / MMsys.MC.barofreq)):
            MMsys.MC_Barostat_step()

    elif simulation_type == "Constant_V":
        MMsys.Poisson_solver_fixed_voltage(Niterations=1, compute_intermediate_forces=True, print_flag=True)
        for j in range(int(freq_traj_output_ps * 1000 / freq_charge_update_fs)):
            MMsys.Poisson_solver_fixed_voltage(Niterations=1, compute_intermediate_forces=True)
            MMsys.simmd.step(freq_charge_update_fs)
        if write_charges:
            with open(charge_name, "a") as chargeFile:
                MMsys.write_electrode_charges(chargeFile)
    else:
        print('simulation type not recognized ...')
        sys.exit()

if simulation_type == "MC_equil":
    state = MMsys.simmd.context.getState(getEnergy=True, getForces=True, getVelocities=False, getPositions=True)
    with open('equilibrated.pdb', 'w') as handle:
        PDBFile.writeFile(MMsys.simmd.topology, state.getPositions(), handle)

print('done!')
sys.exit()
