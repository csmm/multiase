from ase.data import covalent_radii
import numpy as np

# itertools.combinations and itertools.permutations are not in python 2.4

try:
	from itertools import combinations, permutations
except:
	def permutations(iterable, r=None):
		# permutations('ABCD', 2) --> AB AC AD BA BC BD CA CB CD DA DB DC
		# permutations(range(3)) --> 012 021 102 120 201 210
		pool = tuple(iterable)
		n = len(pool)
		if r is None: r = n
		if r > n:
			return
		indices = range(n)
		cycles = range(n, n-r, -1)
		yield tuple(pool[i] for i in indices[:r])
		while n:
			for i in reversed(range(r)):
				cycles[i] -= 1
				if cycles[i] == 0:
					indices[i:] = indices[i+1:] + indices[i:i+1]
					cycles[i] = n - i
				else:
					j = cycles[i]
					indices[i], indices[-j] = indices[-j], indices[i]
					yield tuple(pool[i] for i in indices[:r])
					break
			else:
				return
            
	def combinations(iterable, r):
		pool = tuple(iterable)
		n = len(pool)
		for indices in permutations(range(n), r):
			if sorted(indices) == list(indices):
				yield tuple(pool[i] for i in indices)


class Bonds:
	def __init__(self, atoms, pairs=None, autodetect=False, tolerance=1.15):
		self.len_atoms = len(atoms)
		if autodetect:
			self.detect(atoms, tolerance)
		elif pairs != None:
			self.pairs = np.array(pairs)
		else:
			self.pairs = np.zeros((0,2), dtype=int)
			
	def detect(self, atoms, tolerance):
		frac_coords = atoms.get_scaled_positions()
		covradii = np.array([covalent_radii[n] for n in atoms.numbers])
		
		pairs = []
		for i in range(len(atoms)-1):
		  bondvectors = frac_coords[i+1:] - frac_coords[i]
		  bondvectors -= atoms.pbc * np.round(bondvectors)
		  dist_sq = (np.dot(bondvectors, atoms.cell)**2).sum(axis=1)
		  R1 = covradii[i]
		  R2 = covradii[i+1:]
		  bonded = np.where(dist_sq < ((R1+R2)*tolerance)**2)[0]
		  for j in bonded:
			  pairs.append((i,i+1+j))
			  
		self.pairs = np.array(pairs)
		"""
		distances_sq = np.zeros((len(atoms), len(atoms)))
		for dim in range(3):
			X2, X1 = np.meshgrid(frac_coords[:,dim], frac_coords[:,dim])
			fracdist = X2 - X1
			if atoms.pbc[dim]:
				fracdist -= np.round(fracdist) # MIC
			distances_sq += (fracdist*np.linalg.norm(atoms.cell[dim]))**2
		
		covradii = np.array([covalent_radii[n] for n in atoms.numbers])
		R2, R1 = np.meshgrid(covradii, covradii)
		
		bonded = ((R1+R2)*tolerance)**2 > distances_sq
		bonded = np.triu(bonded, 1)
		self.pairs = np.array(np.where(bonded)).T
		"""
	
	
	def __getitem__(self, i):
		hits = np.where(self.pairs == i)
		if len(hits[0]) == 0: return []
		rows, cols = hits
		cols = 1 - cols
		bonded = self.pairs[rows, cols]
		return list(bonded)
		
	def __iter__(self):
		for pair in self.pairs:
			yield tuple(pair)
			
	def __len__(self):
		return len(self.pairs)
	
	def add(self, i, j):
		self.pairs = np.append(self.pairs, np.array((i,j), ndmin=2), axis=0)
		
	def remove(self, i, j):
		hits = ((self.pairs == (i,j)) + (self.pairs == (j,i)))
		if not np.any(hits):
			raise IndexError('Bond (%i, %i) does not exist!' % (i,j))
		self.pairs = self.pairs[~hits.all(1),:]
	
	def get_bond_matrix(self):
		matrix = np.zeros((self.len_atoms, self.len_atoms))
		matrix[self.pairs[:,0], self.pairs[:,1]] = 1
		return matrix + matrix.T
		
	def atoms_extending(self, new_atoms):
		if 'bonds' in new_atoms.info:
			new_pairs = new_atoms.info['bonds'].pairs + self.len_atoms
			self.pairs = np.append(self.pairs, new_pairs, axis=0)
		self.len_atoms += len(new_atoms)
	
	def find_angles(self):
		angles = []
		for j in range(self.len_atoms):
			for i, k in combinations(self[j], 2):
				angles.append((i,j,k))
		return angles
	
	def find_dihedrals(self):
		dihedrals = []
		sets = list()
		for j, k in self.pairs:
			for i in self[j]:
				if i == k:
					continue
				for l in self[k]:
					if l != j:
						dihedrals.append((i,j,k,l))
		return dihedrals
	
	def find_impropers(self):
		impropers = []
		for i in range(self.len_atoms):
			for j,k,l in combinations(self[i], 3):
				impropers.append((i,j,k,l))
		return impropers
