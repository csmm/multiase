from ase import Atoms
from gpaw import GPAW

d = 0.74
a = 6.0

atoms_all = Atoms("H3",
                positions = [(0, 0, 0),
                (0, 0, d),
                (0, 0, 2*d)],
                cell = (a, a, a))
atoms_sides = Atoms("H2",
                positions = [(0, 0, 0),
                (0, 0, 2*d)],
                cell = (a, a, a))

atoms_left = Atoms("H2",
                positions = [(0, 0, 0),
                (0, 0, d)],
                cell = (a, a, a))

atoms_single = Atoms("H1",
                positions = [(0,0,0)],
                cell = (a, a, a))

atoms_all.center()
atoms_sides.center()
atoms_left.center()
atoms_single.center()

calc_1 = GPAW(nbands=3, txt="h1.txt")
calc_2 = GPAW(nbands=3, txt="h2.txt")
calc_3 = GPAW(nbands=3, txt="h3.txt")
calc_4 = GPAW(nbands=3, txt="h4.txt")

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

