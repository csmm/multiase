from ase import Atoms
from gpaw import GPAW
from gpaw.mpi import rank
from ase.data import s22

from multiasecalc.mixer.selector import AtomListSelector
from multiasecalc.mixer.mixer import Mixer, EnergyCalculation, ForceCalculation
from multiasecalc.lammps.reaxff import ReaxFF
from multiasecalc.utils import get_datafile

import numpy as np


dimer_name = "2-pyridoxine_2-aminopyridine_complex"

atoms = s22.create_s22_system(dimer_name) # full system
atom_counts = s22.get_number_of_dimer_atoms(dimer_name)
Mixer.set_atom_ids(atoms) # sequence of numbers in same order as positions
                          # were given above, index starts from 0

atoms_m1 = atoms[:atom_counts[0]] # first molecule
atoms_m2 = atoms[atom_counts[0]:] # second molecule

calc_gpaw = GPAW(nbands=-2, txt="h2_1.txt")
calc_reaxff = ReaxFF(ff_file_path=get_datafile("ffield.reax.new"),
                implementation="C")
reaxff_cell = (100.0, 100.0, 100.0)
gpaw_cell = (10.0, 10.0, 10.0)

filter_full_system = AtomListSelector(range(len(atoms)),
                                      dict(zip(range(len(atoms)),
                                               [1.0] * len(atoms))))


filter_qm_region = AtomListSelector(range(atom_counts[0]),
                                    dict(zip(range(atom_counts[0]),
                                             [1.0] * atom_counts[0])))

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
mixer = Mixer(name="dimer_mixer", forces=mixer_forces,
              energies=mixer_energies)

atoms.set_calculator(mixer)
atoms_m1.set_calculator(mixer)
atoms_m2.set_calculator(mixer)
m1_energy = atoms_m1.get_potential_energy()
m2_energy = atoms_m2.get_potential_energy()
total_energy = atoms.get_potential_energy()
if rank == 0:
    print("total: %f m1: %f m2: %f interaction: %f" % (
        total_energy, m1_energy, m2_energy,
        total_energy - m1_energy - m2_energy))
