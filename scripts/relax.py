from ase import Atoms, Atom
from ase.io import read
from gpaw import GPAW
from ase.parallel import parprint
from ase.optimize import QuasiNewton


def set_work_cell(molecule,h,vacuum):
    molecule.center(vacuum=vacuum)
    cell = molecule.get_cell()
    for i in [0,1,2]:
        cell[i,i] = int(cell[i,i]/4./h)*4.*h
    molecule.set_cell(cell)
    return

h=0.2
vacuum=4.5
atoms = read('init.traj')
#atoms = read('testMolecule1nm.xyz')
#atoms.set_pbc(True)
#atoms.set_cell([10.0,10.0,10.0])
#set_work_cell(atoms,h,vacuum)

## gpaw calculator:
calc = GPAW(h=h, nbands=-20, xc='PBE', txt='relax.txt',
            spinpol=True, charge=+1)
atoms.set_calculator(calc)

relax = QuasiNewton(atoms, logfile='relax.log',trajectory='relax.traj')
relax.run(fmax=0.05)

calc.write('relax.gpw',mode='all')


