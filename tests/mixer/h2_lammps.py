from ase import Atoms
from csmmcalc.lammps.reaxff import ReaxFF
from csmmcalc.utils import get_datafile

d = 0.74
#d = 3.0
a = 6.0

atoms_all = Atoms("H3",
                positions = [(0, 0, 0),
                (0, 0, d),
                (0, 0, 2*d)],
                cell = (100*a, 100*a, 100*a))
atoms_all.center()

calc_1 = ReaxFF(specorder = ("C", "H", "O", "N", "S"),
        implementation="C",
        ff_file_path=get_datafile("ffield.reax"))

atoms_all.set_calculator(calc_1)
print(atoms_all.get_forces())
print(atoms_all.get_potential_energy())

