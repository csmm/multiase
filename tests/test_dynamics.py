from ase.atoms import Atoms
from multiasecalc.lammps.reaxff import ReaxFF
from multiasecalc.lammps.compass import COMPASS
from multiasecalc.lammps.dynamics import LAMMPSOptimizer, LAMMPS_NVT
from multiasecalc.utils import get_datafile
from ase.data import s22
from ase import units
import numpy as np
import ase.io
from ase.io.trajectory import PickleTrajectory

atoms = s22.create_s22_system('Methane_dimer')
atoms.center(vacuum=10.0)
print atoms.positions

atoms.calc = COMPASS(ff_file_path=get_datafile('compass.frc'))
optimizer = LAMMPSOptimizer(atoms)
optimizer.run()
print atoms.positions

atoms.calc = ReaxFF(ff_file_path=get_datafile('ffield.reax'))
dyn = LAMMPS_NVT(atoms, 1*units.fs, 100, trajectory='test.traj', traj_interval = 2)
dyn.run(5)

atoms.calc = COMPASS(ff_file_path=get_datafile('compass.frc'))
dyn.run(10)

trj = PickleTrajectory('test.traj', 'r')

for t in trj: print t.positions
