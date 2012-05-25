from ase.atoms import Atoms
from lammpsCalculators.reaxff import ReaxFF
from lammpsCalculators.compass import COMPASS
from lammpsCalculators.dynamics import LAMMPSOptimizer, LAMMPS_NVT
from ase.data import s22
from ase import units
import numpy as np
import ase.io

def get_s22_system(name):
	moldata = s22.data[name]
	atoms = Atoms(moldata['symbols'], moldata['positions'])
	
	xmin = min(atoms.positions[:, 0])
	ymin = min(atoms.positions[:, 1])
	zmin = min(atoms.positions[:, 2])
	
	xmax = max(atoms.positions[:, 0])
	ymax = max(atoms.positions[:, 1])
	zmax = max(atoms.positions[:, 2])
	
	atoms.translate(np.array((-xmin, -ymin, -zmin))*2)
	atoms.set_cell(np.array([(xmax-xmin, 0, 0), (0, ymax-ymin, 0), (0, 0, zmax-zmin)])*4)
	return atoms


atoms = get_s22_system('Methane_dimer')
print atoms.positions

atoms.calc = COMPASS()
optimizer = LAMMPSOptimizer(atoms)
optimizer.run()
print atoms.positions

atoms.calc = ReaxFF(specorder = ('C', 'H', 'O', 'N', 'S'))
dyn = LAMMPS_NVT(atoms, 1*units.fs, 100)
dyn.run(1000)
print atoms.positions

atoms.calc = COMPASS()
dyn.run(10000)

ase.io.write('test.xyz', atoms)

print atoms.positions