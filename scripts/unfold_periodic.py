import numpy as np

def unfold_periodic(atoms):
	frac_positions = atoms.get_scaled_positions()
	bonds = atoms.info['bonds']
	
	for i in range(len(atoms)):
		for j in bonds[i]:
			if j > i:
				r = frac_positions[j] - frac_positions[i]
				frac_positions[j] -= np.floor(r + .5)
	
	atoms.set_scaled_positions(frac_positions)
	

if __name__ == '__main__':
	import sys
	from ase.io.trajectory import PickleTrajectory
	from ase.visualize import view
	from multiasecalc.lammps import Bonds
	filename = sys.argv[1]
	traj = PickleTrajectory(filename)
	if len(sys.argv) == 3:
		out = PickleTrajectory(sys.argv[2], 'w')
		for atoms in traj:
			unfold_periodic(atoms)
			out.write(atoms)
	else:
		atoms = traj[-1]
		atoms.info['bonds'] = Bonds(atoms, autodetect=True)
		unfold_periodic(atoms)
		view(atoms)
	