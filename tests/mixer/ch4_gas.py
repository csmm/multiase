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
from csmmcalc.utils import get_datafile
from csmmcalc.mixer.selector import CalcBox

import numpy as np
import cPickle as pickle


# PARAMETERS
a = 20.0
calc_style = "combined" # "combined", "classical", "quantum"
md_style = "Verlet" # "Verlet", "Langevin"
T = 300 # Kelvin
timestep = 0.1*units.fs
# END OF PARAMETERS

# verify that MPI is actually working
print("rank: %i" % rank)



# a cell full of methane

# unpickle atoms, methane_count, cell variables
fd = open("ch4_gas.pkl", "r")
methane_count = pickle.load(fd)
cell = pickle.load(fd)
atoms = pickle.load(fd)
fd.close()

Mixer.set_atom_ids(atoms) # this one is important!

calc_gpaw = GPAW(nbands=-2, txt="h2_1_%i.txt" % rank)
calc_reaxff_full = ReaxFF(ff_file_path=get_datafile("ffield.reax.new"),
                          implementation="C")
calc_reaxff_qbox = ReaxFF(ff_file_path=get_datafile("ffield.reax.new"),
                          implementation="C")

filter_full_sys = CalcBox(pos=(0,0,0), dim=cell, cutoff=2.0, pbc=(1,1,1),
                          debug=0)
filter_qbox = CalcBox(pos=(0,0,0), dim=(a,a,a),
        cutoff=2.0, inner_dim=(a-4,a-4,a-4),
        debug=0)

# full system classical is taken as positive
forces_full_system = ForceCalculation(filter_full_sys)
forces_full_system.calculator = calc_reaxff_full
forces_full_system.cell = cell

# quantum box classical is subtracted using the qbox weights
forces_qbox_reaxff = ForceCalculation(filter_qbox)
forces_qbox_reaxff.calculator = calc_reaxff_qbox
forces_qbox_reaxff.cell = cell
forces_qbox_reaxff.coeff = -1.0

# quantum box quantum is added using qbox weights
forces_qbox_gpaw = ForceCalculation(filter_qbox)
forces_qbox_gpaw.calculator = calc_gpaw
forces_qbox_gpaw.cell = (a+2,a+2,a+2)

# energies are based on H = H_c + H_q' - H_c'
energy_full_system = EnergyCalculation(filter_full_sys)
energy_full_system.calculator = calc_reaxff_full
energy_full_system.cell = cell
energy_full_system.coeff = 1.0

energy_qbox_reaxff = EnergyCalculation(filter_qbox)
energy_qbox_reaxff.calculator = calc_reaxff_qbox
energy_qbox_reaxff.cell = cell
energy_qbox_reaxff.coeff = -1.0

energy_qbox_gpaw = EnergyCalculation(filter_qbox)
energy_qbox_gpaw.calculator = calc_gpaw
energy_qbox_gpaw.cell = (a+2,a+2,a+2)
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

mixer = Mixer(forces=mixer_forces,
              energies=mixer_energies)
atoms.set_calculator(mixer)

dyn = None

if md_style == "Langevin":
    dyn = Langevin(atoms, timestep, 1.5*T*units.kB, 0.002)
elif md_style == "Verlet":
    # set momenta to match temperature T
    MaxwellBoltzmannDistribution(atoms, 1.5*T*units.kB)
    dyn = VelocityVerlet(atoms, timestep)

def printenergy(a=atoms):    #store a reference to atoms in the definition.
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print("Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  Etot = %.3feV" %
        (epot, ekin, ekin/(1.5*units.kB), epot+ekin))

if rank == 0:
    dyn.attach(printenergy, interval=25)
    traj = PickleTrajectory("ch4_gas.traj", 'w', atoms)
    dyn.attach(traj.write, interval=25)
    printenergy()

dyn.run(10000)

