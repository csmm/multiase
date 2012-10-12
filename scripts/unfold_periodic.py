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
	from scripts import build
	from ase.visualize import view
	atoms = build.createPVA(15, 3)
	unfold_periodic(atoms)
	view(atoms)