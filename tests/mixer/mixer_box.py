# Copyright (C) 2012 Aalto University
# Author: Lauri Leukkunen <lauri.leukkunen@aalto.fi>

from ase import Atoms
from ase.constraints import FixBondLengths
from ase.data.molecules import molecule
from ase.md.langevin import Langevin
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.io.trajectory import PickleTrajectory
from ase import units
from gpaw.mpi import rank

import os
import sys

from gpaw import GPAW

from csmmcalc.mixer.mixer import Mixer
from csmmcalc.mixer.mixer import EnergyCalculation, ForceCalculation

from csmmcalc.lammps.reaxff import ReaxFF
from csmmcalc.utils import get_datafile, DynTesting
from csmmcalc.mixer.selector import CalcBox

import numpy as np
import cPickle as pickle

import getopt


def usage():
    print("""
python mixer_box.py [options]
options:
    -h, --help                      help
    -c, --classical                 plain classical simulation
    -q, --quantum                   plain quantum simulation
    -d, --debug=N                   debug level
    -i, --input=input.traj          input trajectory file
    -o, --output=output.traj        output trajectory file
    -B, --box=NxNxN                 full system cell
    -Q, --quantum-box=NxNxN         quantum cell
    -T, --transition-buffer=N       transition buffer in angstroms
    -C, --cutoff=N                  molecule bond detection cutoff in
                                    angstroms
    -S, --steps=N                   number of simulation steps to run
    -t, --time-step=N               time step in femtoseconds
    -M, --molecular-dynamics=name   choose from "Langevin", "Verlet",
                                    "TESTING"
    -P, --position=N,N,N            quantum box center position, the
                                    full system box is always centered
                                    at (0.0, 0.0, 0.0)
    -L, --log-interval=N            write system state to output and
                                    energy log every Nth step
    -H, --langevin-temp=N           temperature target in Kelvins for Langevin
                                    molecular dynamics
    -b, --bands=N                   number of electron bands for GPAW
                                    calculation
    N is float (1.234) for all except --debug, --steps, --log-interval
    and --bands where it is an integer.
    """)

try:
    opts, args = getopt.getopt(sys.argv[1:],
            "H:L:B:Q:T:C:S:t:hcqo:d:i:M:b:",[
                "langevin-temp=",
                "log-interval=",
                "box=",
                "quantum-box=",
                "transition-buffer=",
                "cutoff=",
                "steps=",
                "time-step=",
                "help",
                "classical",
                "quantum",
                "output=",
                "debug=",
                "input=",
                "molecular-dynamics=",
                "bands"])
except getopt.GetoptError, err:
    print(str(err))
    usage()
    sys.exit(2)

# PARAMETERS
s = 10.0
cell = (100.0, 100.0, 100.0)
q_cell = (10.0, 10.0, 10.0)
qbox_pos = (0.0, 0.0, 0.0)
transition_buffer = 2.0
cutoff = 2.0
calc_style = "combined" # "combined", "classical", "quantum"
md_style = "Verlet" # "Verlet", "Langevin"
timestep = 0.1*units.fs
output_file = "mixer_box.traj"
input_file = None
set_debug = 0
steps = 10
log_interval = 25
langevin_temp = 300 # Kelvins
bands = -2
# END OF PARAMETERS

# process command line options
for o, a in opts:
    if o in ["-d", "--debug"]:
        set_debug = int(a)
    if o in ["-h", "--help"]:
        usage()
        sys.exit(0)
    if o in ["-c", "--classical"]:
        calc_style = "classical"
    if o in ["-q", "--quantum"]:
        calc_style = "quantum"
    if o in ["-o", "--output"]:
        output_file = a
    if o in ["-i", "--input"]:
        input_file = a
    if o in ["-M", "--molecular-dynamics"]:
        md_style = a
    if o in ["-t", "--time-step"]:
        timestep = float(a) * units.fs
    if o in ["-B", "--box"]:
        cell = tuple([float(f) for f in a.split("x")])
    if o in ["-Q", "--quantum-box"]:
        q_cell = tuple([float(f) for f in a.split("x")])
    if o in ["-T", "--transition-buffer"]:
        transition_buffer = float(a)
    if o in ["-C", "--cutoff"]:
        cutoff = float(a)
    if o in ["-S", "--steps"]:
        steps = int(a)
    if o in ["-P", "--position"]:
        qbox_pos = tuple([float(f) for f in a.split(",")])
    if o in ["-L", "--log-interval"]:
        log_interval = int(a)
    if o in ["-H", "--langevin-temp"]:
        langevin_temp = float(a)
    if o in ["-b", "--bands"]:
        bands = int(a)


# verify that MPI is actually working
print("rank: %i" % rank)



pt = PickleTrajectory(input_file, "r")
atoms = pt[-1] # get the last step

Mixer.set_atom_ids(atoms) # this one is important!

calc_gpaw = GPAW(nbands=bands, txt="mixer_box_gpaw.log")
calc_reaxff_full = ReaxFF(ff_file_path=get_datafile("ffield.reax.new"),
                          implementation="C")
calc_reaxff_qbox = ReaxFF(ff_file_path=get_datafile("ffield.reax.new"),
                          implementation="C")

# debug disabled for non-master nodes, this is so on purpose!
debug = 0
if rank == 0:
    debug = set_debug

filter_full_sys = CalcBox(name="full_sys",
                pos=(0,0,0), dim=cell,
                cutoff=cutoff, pbc=(1,1,1),
                debug=debug)
filter_qbox = CalcBox(name="qbox",pos=qbox_pos, dim=q_cell,
        cutoff=cutoff,
        inner_dim=(q_cell[0] - transition_buffer*2.0,
            q_cell[1] - transition_buffer*2.0,
            q_cell[2] - transition_buffer*2.0),
        debug=debug)

# full system classical is taken as positive
forces_full_system = ForceCalculation("force_full", filter_full_sys,
                                      debug=debug)
forces_full_system.calculator = calc_reaxff_full
forces_full_system.cell = cell

# quantum box classical is subtracted using the qbox weights
forces_qbox_reaxff = ForceCalculation("force_qbox_reax", filter_qbox,
                                      debug=debug)
forces_qbox_reaxff.calculator = calc_reaxff_qbox
forces_qbox_reaxff.cell = cell
forces_qbox_reaxff.coeff = -1.0

# quantum box quantum is added using qbox weights
forces_qbox_gpaw = ForceCalculation("force_qbox_gpaw", filter_qbox,
                                    debug=debug)
forces_qbox_gpaw.calculator = calc_gpaw
forces_qbox_gpaw.cell = (q_cell[0] + 3,
                         q_cell[1] + 3,
                         q_cell[2] + 3)

# energies are based on H = H_c + H_q' - H_c'
energy_full_system = EnergyCalculation("energy_full", filter_full_sys)
energy_full_system.calculator = calc_reaxff_full
energy_full_system.cell = cell
energy_full_system.coeff = 1.0

energy_qbox_reaxff = EnergyCalculation("energy_qbox_reax", filter_qbox)
energy_qbox_reaxff.calculator = calc_reaxff_qbox
energy_qbox_reaxff.cell = cell
energy_qbox_reaxff.coeff = -1.0

energy_qbox_gpaw = EnergyCalculation("energy_qbox_gpaw", filter_qbox)
energy_qbox_gpaw.calculator = calc_gpaw
energy_qbox_gpaw.cell = (q_cell[0] + 3.0,
                         q_cell[1] + 3.0,
                         q_cell[2] + 3.0)
energy_qbox_gpaw.coeff = 1.0

mixer_forces = []
mixer_energies = []

if calc_style == "combined":
    mixer_forces = [forces_full_system,
                    forces_qbox_reaxff,
                    forces_qbox_gpaw]

    mixer_energies = [energy_full_system,
                      energy_qbox_reaxff,
                      energy_qbox_gpaw]
elif calc_style == "classical":
    mixer_forces = [forces_full_system]
    mixer_energies = [energy_full_system]
elif calc_style == "quantum":
    mixer_forces = [forces_qbox_gpaw]
    mixer_energies = [energy_qbox_gpaw]

mixer = Mixer(name="mixer_box",
              forces=mixer_forces,
              energies=mixer_energies,
              debug=debug)
atoms.set_calculator(mixer)

dyn = None

if md_style == "Langevin":
    dyn = Langevin(atoms, timestep, 1.5*T*units.kB, 0.002)
elif md_style == "Verlet":
    dyn = VelocityVerlet(atoms, timestep)
elif md_style == "TESTING":
    dyn = DynTesting(atoms, timestep=0.1*units.fs, offset=(0.1, 0., 0.))

energy_file = None

def printenergy(a=atoms):
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    energy_file.write("%.5e,%.5e,%.5e,%.5e" %
        (epot + ekin, epot, ekin, ekin/(1.5*units.kB)) + "\n")

# only enable logging for master node where rank == 0
if rank == 0:
    energy_file = open("mixer_box_energy.csv", "w+b", buffering=0)
    dyn.attach(printenergy, interval=log_interval)
    traj = PickleTrajectory(output_file, 'w', atoms)
    dyn.attach(traj.write, interval=log_interval)
    printenergy()

dyn.run(steps)

