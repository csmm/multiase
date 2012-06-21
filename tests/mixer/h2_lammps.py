from ase import Atoms
from csmmcalc.lammps.reaxff import ReaxFF
from csmmcalc.utils import get_datafile

d = 0.74
a = 6.0

atoms = Atoms("H2",
                positions = [(0, 0, 0),
                (0, 0, d)],
                cell = (100*a, 100*a, 100*a))

atoms.center()
calc = ReaxFF(specorder = ("C", "H", "O", "N", "S"),
        ff_file_path=get_datafile("ffield.reax"))
calc.keep_tmp_files = True
atoms.set_calculator(calc)
print(atoms.get_forces())

