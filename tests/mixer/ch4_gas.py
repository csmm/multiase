# Copyright (C) 2012 Aalto University
# Author: Lauri Leukkunen <lauri.leukkunen@aalto.fi>

from ase import Atoms
from ase.constraints import FixBondLengths
from ase.data.molecules import molecule
from ase.md.langevin import Langevin
from ase.io.trajectory import PickleTrajectory
from ase import units

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

d = 0.76470
#d = 1.5
a = 15.0



# a cell full of methane randomly positioned

# unpickle atoms, methane_count, cell variables
fd = open("ch4_gas.pkl", "r")
methane_count = pickle.load(fd)
cell = pickle.load(fd)
atoms = pickle.load(fd)
fd.close()

Mixer.set_atom_ids(atoms) # this one is important!

calc_gpaw = GPAW(nbands=-2, txt="h2_1.txt")
calc_reaxff_full = ReaxFF(ff_file_path=get_datafile("ffield.reax.new"))
calc_reaxff_qbox = ReaxFF(ff_file_path=get_datafile("ffield.reax.new"))

filter_full_sys = CalcBox(pos=(0,0,0), dim=cell, pbc=(1,1,1))
filter_qbox = CalcBox(pos=(0,0,0), dim=(a,a,a), inner_dim=(a-4,a-4,a-4))

# full system classical is taken as positive
forces_full_system = ForceCalculation(filter_full_sys)
forces_full_system.calculator = calc_reaxff_full
forces_full_system.cell = cell

# quantum box classical is deducted using the qbox weights
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


# doing just the classical for now

mixer_forces = [forces_full_system,
                forces_qbox_reaxff,
                forces_qbox_gpaw]
#mixer_forces = [forces_full_system]

mixer_energies = [energy_full_system,
                  energy_qbox_reaxff,
                  energy_qbox_gpaw]
#mixer_energies = [energy_full_system]


mixer = Mixer(forces=mixer_forces,
              energies=mixer_energies)
atoms.set_calculator(mixer)

T = 300 # Kelvin

dyn = Langevin(atoms, 0.1*units.fs, T*units.kB, 0.002)

def printenergy(a=atoms):    #store a reference to atoms in the definition.
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print ("Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  Etot = %.3feV" %
        (epot, ekin, ekin/(1.5*units.kB), epot+ekin))

dyn.attach(printenergy, interval=5)
traj = PickleTrajectory("ch4_gas.traj", 'w', atoms)
dyn.attach(traj.write, interval=5)
printenergy()
dyn.run(10000)

