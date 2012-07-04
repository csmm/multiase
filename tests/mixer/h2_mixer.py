from ase import Atoms
from gpaw import GPAW

from csmmcalc.mixer.mixer import Mixer
from csmmcalc.mixer.mixer import EnergyCalculation, ForceCalculation
from csmmcalc.lammps.reaxff import ReaxFF
from csmmcalc.utils import get_datafile

import numpy as np

d = 0.74
#d = 1.5
a = 6.0

atoms = Atoms("H3",
                positions = [(0, 0, 0),
                (0, 0, d),
                (0, 0, 2*d)],
                cell = (10*a, 10*a, 10*a))

Mixer.set_atom_ids(atoms) # sequence of numbers in same order as positions
                          # were given above, index starts from 0

calc_1 = GPAW(nbands=3, txt="h2_1.txt")
#calc_2 = GPAW(nbands=3, txt="h2_2.txt")
calc_3 = ReaxFF(ff_file_path=get_datafile("ffield.reax.new"))


# the keys below must match tags given above
forces_full_system = ForceCalculation()
forces_full_system.calculator = calc_3
forces_full_system.atom_ids = Mixer.get_atom_ids(atoms)
forces_full_system.cell = (100*a, 100*a, 100*a)
forces_full_system.weights = {0: 0.25,
                              1: 0.5,
                              2: 0.75}

forces_qm_region = ForceCalculation()
forces_qm_region.calculator = calc_1
forces_qm_region.atom_ids = (0, 1, 2)
forces_qm_region.cell = (a, a, a)
forces_qm_region.weights = {0: 0.75,
                            1: 0.5,
                            2: 0.25}

energy_full_system_reaxff = EnergyCalculation()
energy_full_system_reaxff.calculator = calc_3
energy_full_system_reaxff.atom_ids = Mixer.get_atom_ids(atoms)
energy_full_system_reaxff.cell = (100*a, 100*a, 100*a)
energy_full_system_reaxff.coeff = 1

energy_qm_region_reaxff = EnergyCalculation()
energy_qm_region_reaxff.calculator = calc_3
energy_qm_region_reaxff.atom_ids = (0, 1)
energy_qm_region_reaxff.cell = (100*a, 100*a, 100*a)
energy_qm_region_reaxff.coeff = -1

energy_qm_region_gpaw = EnergyCalculation()
energy_qm_region_gpaw.calculator = calc_1
energy_qm_region_gpaw.atom_ids = (0, 1)
energy_qm_region_gpaw.cell = (a, a, a)
energy_qm_region_gpaw.coeff = 1

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
