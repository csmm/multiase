from ase import Atoms
from ase.units import Bohr, Hartree
from gpaw import GPAW
from ase.parallel import rank
from ase.parallel import parprint

d = 1.4 * Bohr
#print("d: %f" % d)
a = 6.0

atoms_all = Atoms("H2",
                positions = [(0, 0, 0),
                (0, 0, d)],
                cell = (a, a, a))

atoms_all.center()

calc_1 = GPAW(nbands=2, txt="h1.txt")

atoms_all.set_calculator(calc_1)

parprint(atoms_all.get_forces())
parprint(atoms_all.get_potential_energy())
parprint("H2 energy in eV: %f" % (-1.174475*Hartree))
