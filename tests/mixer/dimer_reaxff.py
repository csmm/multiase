from ase import Atoms
from gpaw import GPAW
from ase.data import s22

from multiasecalc.mixer.selector import AtomListSelector
from multiasecalc.mixer.mixer import Mixer, EnergyCalculation, ForceCalculation
from multiasecalc.lammps.reaxff import ReaxFF
from multiasecalc.utils import get_datafile

import numpy as np


dimer_name = "2-pyridoxine_2-aminopyridine_complex"

atoms = s22.create_s22_system(dimer_name)
atom_counts = s22.get_number_of_dimer_atoms(dimer_name)
atoms_m1 = atoms[:atom_counts[0]] # first molecule
atoms_m2 = atoms[atom_counts[0]:] # second molecule

calc_reaxff = ReaxFF(ff_file_path=get_datafile("ffield.reax.new"),
                implementation="C")
reaxff_cell = (100.0, 100.0, 100.0)
atoms.set_cell(reaxff_cell)
atoms_m1.set_cell(reaxff_cell)
atoms_m2.set_cell(reaxff_cell)
atoms.center()
atoms_m1.center()
atoms_m2.center()
atoms_m1.set_calculator(calc_reaxff)
atoms_m2.set_calculator(calc_reaxff)
atoms.set_calculator(calc_reaxff)

m1_energy = atoms_m1.get_potential_energy()
m2_energy = atoms_m2.get_potential_energy()
total_energy = atoms.get_potential_energy()
print("total: %f m1: %f m2: %f interaction: %f" % (
        total_energy, m1_energy, m2_energy,
        total_energy - m1_energy - m2_energy))
