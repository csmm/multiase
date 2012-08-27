from ase import Atoms
from gpaw import GPAW
from ase.data import s22
from gpaw.mpi import rank

from csmmcalc.mixer.selector import AtomListSelector
from csmmcalc.mixer.mixer import Mixer, EnergyCalculation, ForceCalculation
from csmmcalc.lammps.reaxff import ReaxFF
from csmmcalc.utils import get_datafile

import numpy as np


dimer_name = "2-pyridoxine_2-aminopyridine_complex"

atoms = s22.create_s22_system(dimer_name)
atom_counts = s22.get_number_of_dimer_atoms(dimer_name)
atoms_m1 = atoms[:atom_counts[0]] # first molecule
atoms_m2 = atoms[atom_counts[0]:] # second molecule


calc_gpaw = GPAW(nbands=-2, xc="PBE", txt="dimer_gpaw.gpaw.log")
gpaw_cell = (15.0, 15.0, 15.0)
atoms.set_cell(gpaw_cell)
atoms_m1.set_cell(gpaw_cell)
atoms_m2.set_cell(gpaw_cell)
atoms.center()
atoms_m1.center()
atoms_m2.center()
atoms.set_calculator(calc_gpaw)
atoms_m1.set_calculator(calc_gpaw)
atoms_m2.set_calculator(calc_gpaw)

m1_energy = atoms_m1.get_potential_energy()
m2_energy = atoms_m2.get_potential_energy()
total_energy = atoms.get_potential_energy()
if rank == 0:
    print("total: %f m1: %f m2: %f interaction: %f" % (
        total_energy, m1_energy, m2_energy,
        total_energy - m1_energy - m2_energy))
