from ase.data import covalent_radii
import numpy as np

# itertools.combinations and itertools.permutations are not in python 2.4

# from itertools import combinations, permutations

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
	def __init__(self, atoms, pairs=None, autodetect=False):
		self.len_atoms = len(atoms)
		if autodetect:
			self.detect(atoms)
		elif pairs:
			self.pairs = np.array(pairs)
		else:
			self.pairs = np.zeros((0,2), dtype=int)
			
	def detect(self, atoms):
		tolerance = 1.1
		def bonded(i, j):
			bondLength = covalent_radii[atoms[i].number] + covalent_radii[atoms[j].number]
			return atoms.get_distance(i, j, mic=True) < bondLength*tolerance
		
		N = len(atoms)
		pairs = [pair for pair in combinations(range(N), 2) if bonded(*pair)]
		self.pairs = np.array(pairs)
	
	
	def __getitem__(self, i):
		rows, cols = np.where(self.pairs == i)
		cols = 1 - cols
		bonded = self.pairs[rows, cols]
		return tuple(bonded)
		
	def __iter__(self):
		for pair in self.pairs:
			yield tuple(pair)
			
	def __len__(self):
		return len(self.pairs)
		
	def add(self, i, j):
		self.pairs = np.append(self.pairs, np.array((i,j), ndmin=2), axis=0)
		
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
		for j in range(self.len_atoms):
			for i, k in permutations(self[j], 2):
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
