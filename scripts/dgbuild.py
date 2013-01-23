from __future__ import division

from ase import Atoms, Atom
from multiasecalc.lammps import Bonds
import numpy as np
from numpy import pi, cos, sin, sqrt, dot
import random
random.seed(None)
np.random.seed(None)
from itertools import combinations, product

def distance_geometry(mol, bond_matrix, unit_indices, characteristic_ratio=8.3):
	N = len(mol)
	
	BOND_LENGTH = 1.6
	BOND_ANGLE = 113*np.pi/180
	DIST_1_3 = BOND_LENGTH * sqrt(2 - 2*cos(BOND_ANGLE))
	MIN_DIST_1_4 = BOND_LENGTH * (1 + 2*cos(pi-BOND_ANGLE))
	MAX_DIST_1_4 = 2 * sqrt(BOND_LENGTH**2 + (BOND_LENGTH/2)**2 - 2 * BOND_LENGTH**2 /2 * cos(BOND_ANGLE))
	MIN_VDW_DIST = 1.6*2
	MAX_DIST = 99
	
	print 'Create distance matrix'
	
	"""bond_matrix = np.zeros((N,N), dtype = bool)
	for i in range(N):
		for j in bonds[i]:
			bond_matrix[i,j] = True
	"""
	bonds12 = bond_matrix
	bonds13 = dot(bonds12, bond_matrix)
	bonds14 = dot(bonds13, bond_matrix)
		
	Dmax = np.ones((N, N))*MAX_DIST
	Dmin = np.ones((N, N))*MIN_VDW_DIST
	
	Dmax[np.where(bonds14)] = MAX_DIST_1_4
	Dmin[np.where(bonds14)] = MIN_DIST_1_4
	Dmax[np.where(bonds13)] = DIST_1_3
	Dmin[np.where(bonds13)] = DIST_1_3
	Dmax[np.where(bonds12)] = BOND_LENGTH
	Dmin[np.where(bonds12)] = BOND_LENGTH
	
	D = np.triu(Dmax, k=1) + np.tril(Dmin, k=-1)
	
	print D
	
	print 'Triangle smoothing'
	triangle_smooth(D)
	#print D
	
	print 'Generate coarse chain'
	N_monomers = len(mol[unit_indices])
	kuhn_length = BOND_LENGTH * characteristic_ratio / cos(BOND_ANGLE/2)
	N_coarse = N_monomers * cos(BOND_ANGLE/2)**2 / characteristic_ratio
	coarse_segment_units = int(N_monomers/N_coarse)
	print 'Coarse chain length:', N_coarse
	print 'kuhn length =', kuhn_length
	print coarse_segment_units, 'polymer units in a coarse unit'
	
	coarse_indices = slice(unit_indices.start, None, unit_indices.step*coarse_segment_units)
	D_coarse = D[coarse_indices, coarse_indices]
	r_mean = kuhn_length
	gaussian_matrix = gaussian_chain(D_coarse, r_mean)
	print gaussian_matrix
	D[coarse_indices, coarse_indices] = gaussian_matrix
	
	triangle_smooth(D)
	
	print 'Random generation'
	for i, j in combinations(range(N), 2):
		val = random.uniform(D[i,j], D[j,i])
		D[i,j] = val
		D[j,i] = val
	#print D
	
	print 'Metric matrix'
	d0 = D[0,:]
	#print np.mean(D**2, 1) - 1.0/N**2 * np.sum(np.triu(D, k=1)**2)
	#d0 = sqrt( np.mean(D**2, 1) - 1.0/N**2 * np.sum(np.triu(D, k=1)**2) )
	
	d0h, d0v = np.meshgrid(d0, d0)
	G = (d0h**2 + d0v**2 - D**2) / 2
	#print G
	
	print 'Eigenvalues of metric matrix'
	evals, evecs = np.linalg.eig(G)
	
	descending_order = np.argsort(evals)[::-1]
	evals = evals[descending_order]
	evecs = evecs[:,descending_order]
	#print evals
	#print 'Eigenvectors'
	#print evecs
	
	evals3 = np.append(evals[:3], np.zeros(len(evals)-3))
	X = dot(evecs[:,:3], np.diag(sqrt(evals[:3])))
	mol.positions = X
	

def triangle_smooth(D):
	N = D.shape[0]
	Dmax = np.triu(D)
	Dmax = Dmax + Dmax.T
	Dmin = np.tril(D)
	Dmin = Dmin + Dmin.T
	
	changed = True
	
	
	while changed:
		#print 'upper interation'
		changed = False
		
		for i in range(1,N):
			row = Dmax[i,:]
			M = Dmax + np.repeat(np.reshape(row, (N,1)), N, axis=1)
			new_row = np.min(M, axis=0)
			if any(new_row < row):
				new_row = np.min([row, new_row], axis=0)
				changed = True
				Dmax[i,:] = new_row
				Dmax[:,i] = new_row
	
	"""
	while changed:
		changed = False
		for a, b in combinations(range(N), 2):
			bound = np.min(Dmax[a,:]+Dmax[b,:])
			if bound < Dmax[a,b]:
				Dmax[a,b] = bound
				changed = True
	"""
	"""
		for a, b, c in combinations(range(N), 3):
			acmax = D[a,b] + D[b,c]
			if D[a,c] > acmax:
				D[a,c] = acmax
				changed = True
			
			abmax = D[a,c] + D[b,c]
			if D[a,b] > abmax:
				D[a,b] = abmax
				changed = True
				
			bcmax = D[a,b] + D[a,c]
			if D[b,c] > bcmax:
				D[b,c] = bcmax
				changed = True
		"""
		
	changed = True
	
	
	while changed:
		#print 'lower iteration'
		changed = False
		for i in range(1,N):
			maxrow = Dmax[i,:]
			minrow = Dmin[i,:]
			
			M1 = Dmin - np.repeat(np.reshape(maxrow, (N,1)), N, axis=1)
			M2 = np.repeat(np.reshape(minrow, (N,1)), N, axis=1) - Dmax
			new_row1 = np.max(M1, axis=0)
			new_row2 = np.max(M2, axis=0)
			new_row = np.max([new_row1, new_row2, minrow], axis=0)
			if any(new_row != minrow):
				changed = True
				Dmin[i,:] = new_row
				Dmin[:,i] = new_row
	
	"""
	while changed:
		changed = False
		for a, b in combinations(range(N), 2):
			bound = np.max([Dmin[a,:]-Dmax[b,:], Dmin[b,:]-Dmax[a,:]])
			if bound > Dmin[b,a]:
				Dmin[b,a] = bound
				changed = True
	"""
	"""
	while changed:
		changed = False
		for a, b, c in combinations(range(N), 3):
			acmin = max(D[b,a] - D[b,c], D[c,b] - D[a,b])
			if D[c,a] < acmin:
				D[c,a] = acmin
				changed = True
				
			abmin = max(D[c,a] - D[b,c], D[c,b] - D[a,c])
			if D[b,a] < abmin:
				D[b,a] = abmin
				changed = True
				
			bcmin = max(D[b,a] - D[a,c], D[c,a] - D[a,b])
			if D[c,b] < bcmin:
				D[c,b] = bcmin
				changed = True
	"""	
	D[:] = np.tril(Dmin) + np.triu(Dmax)

def gaussian_chain(dist_matrix, r_mean):
	N = dist_matrix.shape[0]
	pos = np.zeros((3))
	points = np.zeros((N, 3))
	for i in range(1,N):
		ok = False
		while not ok:
			step = np.random.normal(0, r_mean/sqrt(3), size=3)
			points[i,:] = points[i-1,:] + step
			D_new = distance_matrix(points)[:i+1,:i+1]
			too_far = np.any(np.triu(D_new) > np.triu(dist_matrix[:i+1,:i+1]))
			too_close = np.any(np.tril(D_new) < np.tril(dist_matrix[:i+1,:i+1]))
			ok = not too_far and not too_close
	return D_new

def distance_matrix(coordinates):
	metric = dot(coordinates, coordinates.T)
	sq_pos = np.diag(metric)
	sq_pos_x, sq_pos_y = np.meshgrid(sq_pos, sq_pos)
	D2 = sq_pos_x + sq_pos_y - 2*metric
	return sqrt(D2)
	
def create_PVA(units):
	start_group = Atoms('H')
	start_bonds = np.matrix('1 1; 1 1')
	
	unit = Atoms('CH2CHOH')
	unit_bonds = np.matrix("""
		1  1  1  1  0  0  0  0 ;
		1  1  0  0  0  0  0  0 ;
		1  0  1  0  0  0  0  0 ;
		1  0  0  1  1  1  0  1 ;
		0  0  0  1  1  0  0  0 ;
		0  0  0  1  0  1  1  0 ;
		0  0  0  0  0  1  1  0 ;
		0  0  0  1  0  0  0  0 
		""")
	unit_carbons = np.array([0, 3])
		
	end_group = Atoms('CH3')
	end_bonds = np.matrix("""
		1 1 1 1 ;
		1 1 0 0 ;
		1 0 1 0 ;
		1 0 0 1 
		""")
	
	N = len(start_group) + units*len(unit) + len(end_group)
	atoms = start_group.copy()
	bonds = np.zeros((N,N), dtype=int)
	bonds[:2,:2] = start_bonds
	
	for i in range(units):
		n = len(atoms)
		atoms.extend(unit)
		bonds[n:n+len(unit)+1, n:n+len(unit)+1] = unit_bonds
	
	n = len(atoms)
	atoms.extend(end_group)
	bonds[n:n+len(end_group), n:n+len(end_group)] = end_bonds
	
	unit_indices = slice(len(start_group), None, len(unit))

	return atoms, bonds, unit_indices

def build_PVA(N):
	mol, bond_matrix, backbone_carbons = create_PVA(N)
	distance_geometry(mol, bond_matrix, backbone_carbons)

	N = len(mol)
	bonds = np.array(np.where(np.triu(bond_matrix, 1))).T
	mol.info['bonds'] = Bonds(mol, pairs=bonds)
	
	from multiasecalc.lammps.compass import COMPASS
	from multiasecalc.lammps.dynamics import LAMMPSOptimizer
	from multiasecalc.utils import get_datafile
	mol.calc = COMPASS(get_datafile('compass.frc'), parameters=dict(extra_cmds=['communicate single cutoff 80']), debug=True)
	dyn = LAMMPSOptimizer(mol)
	dyn.run()
	
	return mol

if __name__ == '__main__':	
	from ase.visualize import view
	mol = build_PVA(20)
	view(mol)