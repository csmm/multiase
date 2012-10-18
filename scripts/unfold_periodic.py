import numpy as np

def unfold_periodic(atoms):
	frac_positions = atoms.get_scaled_positions()
	bonds = atoms.info['bonds']
	
	for i in range(len(atoms)):
		for j in bonds[i]:
			if j > i:
				r = frac_positions[j] - frac_positions[i]
				frac_positions[j] -= np.round(r)
	
	atoms.set_scaled_positions(frac_positions)
	

if __name__ == '__main__':
	import sys
	from ase import io
	from ase.io.trajectory import PickleTrajectory
	from ase.visualize import view
	from multiasecalc.lammps import Bonds
	filename = sys.argv[1]
	if len(sys.argv) == 3:
		traj = PickleTrajectory(filename)
		out = PickleTrajectory(sys.argv[2], 'w')
		for atoms in traj:
			unfold_periodic(atoms)
			out.write(atoms)
	else:
		atoms = io.read(filename)
		if not 'bonds' in atoms.info:
			atoms.info['bonds'] = Bonds(atoms, autodetect=True)
		unfold_periodic(atoms)
		view(atoms)
	