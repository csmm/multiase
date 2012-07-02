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
atoms_sides = Atoms("H2",
                positions = [(0, 0, 0),
                (0, 0, 2*d)],
                cell = (100*a, 100*a, 100*a))

atoms_left = Atoms("H2",
                positions = [(0, 0, 0),
                (0, 0, d)],
                cell = (100*a, 100*a, 100*a))

atoms_single = Atoms("H1",
                positions = [(0,0,0)],
                cell = (100*a, 100*a, 100*a))

atoms_all.center()
atoms_sides.center()
atoms_left.center()
atoms_single.center()

calc_1 = ReaxFF(specorder = ("C", "H", "O", "N", "S"),
        ff_file_path=get_datafile("ffield.reax"))
calc_2 = ReaxFF(specorder = ("C", "H", "O", "N", "S"),
        ff_file_path=get_datafile("ffield.reax"))
calc_3 = ReaxFF(specorder = ("C", "H", "O", "N", "S"),
        ff_file_path=get_datafile("ffield.reax"))
calc_4 = ReaxFF(specorder = ("C", "H", "O", "N", "S"),
        ff_file_path=get_datafile("ffield.reax"))

atoms_all.set_calculator(calc_1)
atoms_sides.set_calculator(calc_2)
atoms_left.set_calculator(calc_3)
atoms_single.set_calculator(calc_4)

print(atoms_all.get_forces())
print(atoms_all.get_potential_energy())

print(atoms_sides.get_forces())
print(atoms_sides.get_potential_energy())

print(atoms_left.get_forces())
print(atoms_left.get_potential_energy())

print(atoms_single.get_forces())
print(atoms_single.get_potential_energy())

