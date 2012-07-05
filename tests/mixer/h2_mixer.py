from ase import Atoms
from gpaw import GPAW

from csmmcalc.mixer.mixer import Mixer
from csmmcalc.mixer.mixer import EnergyCalculation, ForceCalculation
from csmmcalc.lammps.reaxff import ReaxFF
from csmmcalc.utils import get_datafile

import numpy as np

d = 0.76470
#d = 1.5
a = 6.0

atoms = Atoms("H2",
                positions = [(0, 0, 0),
                (0, 0, d)],
                cell = (10*a, 10*a, 10*a))

Mixer.set_atom_ids(atoms) # sequence of numbers in same order as positions
                          # were given above, index starts from 0

calc_1 = GPAW(nbands=2, txt="h2_1.txt")
#calc_2 = GPAW(nbands=2, txt="h2_2.txt")
#calc_2_cell = (a, a, a)
calc_2 = ReaxFF(ff_file_path=get_datafile("ffield.reax.new"))
calc_2_cell = (100*a, 100*a, 100*a)

forces_full_system = ForceCalculation()
forces_full_system.calculator = calc_2
forces_full_system.atom_ids = Mixer.get_atom_ids(atoms)
forces_full_system.cell = calc_2_cell
forces_full_system.weights = {0: 0.0,
                              1: 1.0}

forces_qm_region = ForceCalculation()
forces_qm_region.calculator = calc_1
forces_qm_region.atom_ids = [0, 1]
forces_qm_region.cell = (a, a, a)
forces_qm_region.weights = {0: 1.0,
                            1: 0.0}

energy_full_system_reaxff = EnergyCalculation()
energy_full_system_reaxff.calculator = calc_2
energy_full_system_reaxff.atom_ids = Mixer.get_atom_ids(atoms)
energy_full_system_reaxff.cell = calc_2_cell
energy_full_system_reaxff.coeff = 1.0

energy_qm_region_reaxff = EnergyCalculation()
energy_qm_region_reaxff.calculator = calc_2
energy_qm_region_reaxff.atom_ids = [0]
energy_qm_region_reaxff.cell = calc_2_cell
energy_qm_region_reaxff.coeff = -1.0

energy_qm_region_gpaw = EnergyCalculation()
energy_qm_region_gpaw.calculator = calc_1
energy_qm_region_gpaw.atom_ids = [0]
energy_qm_region_gpaw.cell = (a, a, a)
energy_qm_region_gpaw.coeff = 1.0

mixer_forces = (forces_full_system, forces_qm_region)
mixer_energies = (energy_full_system_reaxff,
                energy_qm_region_reaxff,
                energy_qm_region_gpaw)


atoms.center()
mixer = Mixer(forces=mixer_forces,
              energies=mixer_energies)

atoms.set_calculator(mixer)
print(atoms.get_forces())
print(atoms.get_potential_energy())
