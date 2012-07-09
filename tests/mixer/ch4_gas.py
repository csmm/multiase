# Copyright (C) 2012 Aalto University
# Author: Lauri Leukkunen <lauri.leukkunen@aalto.fi>

from ase import Atoms
from ase.constraints import FixBondLengths
from ase.data.molecules import molecule
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
a = 6.0



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

filter_full_sys = CalcBox(position=(0,0,0), dim=cell, pbc=(1,1,1))
filter_qbox = CalcBox(position=(0,0,0), dim=(a,a,a), inner_dim=(a-2,a-2,a-2))

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

mixer_forces = [forces_full_system,
                forces_qbox_reaxff,
                forces_qbox_gpaw]
mixer_energies = [energy_full_system,
                  energy_qbox_reaxff,
                  energy_qbox_gpaw]


mixer = Mixer(forces=mixer_forces,
              energies=mixer_energies)
atoms.set_calculator(mixer)
print(atoms.get_forces())
print(atoms.get_potential_energy())

