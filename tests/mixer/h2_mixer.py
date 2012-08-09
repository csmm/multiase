from ase import Atoms
from gpaw import GPAW

from csmmcalc.mixer.selector import AtomListSelector
from csmmcalc.mixer.mixer import Mixer, EnergyCalculation, ForceCalculation
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
calc_2 = ReaxFF(ff_file_path=get_datafile("ffield.reax.new"),
                implementation="C")
calc_2_cell = (100*a, 100*a, 100*a)


filter_full_system = AtomListSelector(Mixer.get_atom_ids(atoms),
                        {0: 0.0,
                         1: 1.0})

filter_qm_region = AtomListSelector((0, 1), {0: 1.0, 1: 0.0})

forces_full_system = ForceCalculation("forces_full_sys", filter_full_system)
forces_full_system.calculator = calc_2
forces_full_system.cell = calc_2_cell

forces_qm_region = ForceCalculation("forces_qm", filter_qm_region)
forces_qm_region.calculator = calc_1
forces_qm_region.cell = (a, a, a)

energy_full_system_reaxff = EnergyCalculation("energy_full_sys",
                                              filter_full_system)
energy_full_system_reaxff.calculator = calc_2
energy_full_system_reaxff.cell = calc_2_cell
energy_full_system_reaxff.coeff = 1.0

energy_qm_region_reaxff = EnergyCalculation("energy_qm_reaxff",
                                            filter_qm_region)
energy_qm_region_reaxff.calculator = calc_2
energy_qm_region_reaxff.cell = calc_2_cell
energy_qm_region_reaxff.coeff = -1.0

energy_qm_region_gpaw = EnergyCalculation("energy_qm_gpaw",
                                          filter_qm_region)
energy_qm_region_gpaw.calculator = calc_1
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
