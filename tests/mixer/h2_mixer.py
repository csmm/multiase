from ase import Atoms
from gpaw import GPAW

from multiasecalc.mixer.selector import AtomListSelector
from multiasecalc.mixer.mixer import Mixer, EnergyCalculation, ForceCalculation
from multiasecalc.lammps.reaxff import ReaxFF
from multiasecalc.utils import get_datafile

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

calc_gpaw = GPAW(nbands=2, txt="h2_1.txt")
calc_reaxff = ReaxFF(ff_file_path=get_datafile("ffield.reax.new"),
                implementation="C")
reaxff_cell = (100*a, 100*a, 100*a)
gpaw_cell = (a, a, a)


filter_full_system = AtomListSelector((0, 1),
                        {0: 1.0,
                         1: 1.0})

filter_qm_region = AtomListSelector((0, 1), {0: 1.0, 1: 0.0})

forces_full_system = ForceCalculation("forces_full_sys",
                                      selector=filter_full_system,
                                      calculator=calc_reaxff,
                                      cell=reaxff_cell)

forces_qm_gpaw = ForceCalculation("forces_qm_gpaw",
                                  selector=filter_qm_region,
                                  calculator=calc_gpaw,
                                  cell=gpaw_cell)

forces_qm_reaxff = ForceCalculation("forces_qm_reaxff",
                                    selector=filter_qm_region,
                                    calculator=calc_reaxff,
                                    cell=reaxff_cell)

energy_full_system_reaxff = EnergyCalculation("energy_full_sys",
                                              selector=filter_full_system,
                                              calculator=calc_reaxff,
                                              cell=reaxff_cell)

energy_qm_region_reaxff = EnergyCalculation("energy_qm_reaxff",
                                            selector=filter_qm_region,
                                            calculator=calc_reaxff,
                                            cell=reaxff_cell,
                                            coeff=-1.0)

energy_qm_region_gpaw = EnergyCalculation("energy_qm_gpaw",
                                          selector=filter_qm_region,
                                          calculator=calc_gpaw,
                                          cell=gpaw_cell)

mixer_forces = (forces_full_system, forces_qm_gpaw, forces_qm_reaxff)
mixer_energies = (energy_full_system_reaxff,
                energy_qm_region_reaxff,
                energy_qm_region_gpaw)


atoms.center()
mixer = Mixer(name="H2_mixer", forces=mixer_forces,
              energies=mixer_energies)

atoms.set_calculator(mixer)
print(atoms.get_forces())
print(atoms.get_potential_energy())
