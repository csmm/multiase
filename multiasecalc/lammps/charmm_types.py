from ase import Atoms

funcDict = {}

class CustomAtom(object):
	def __init__(self, symbol):
		self.symbol = symbol
		self.bonded_atoms = []
		self.environment = []
		self.rings = []
		
	def add_bonded(self, atom):
		self.bonded_atoms.append(atom)
		self.environment.append(atom.symbol)
	def aromatic(self):
		return self.in_6_mem_ring() or self.in_5_mem_ring()
	def in_6_mem_ring(self):
		return any([len(ring) == 6 for ring in self.rings])
	def in_5_mem_ring(self):
		return any([len(ring) == 5 for ring in self.rings])
	def heteroatom(self):
		return not self.symbol in ('C', 'H')
	def neighbors_of_type(self, symbol):
		return [n for n in self.bonded_atoms if n.symbol == symbol]
	def n_neigbors(self):
		return len(self.bonded_atoms)
	def has_environment(self, *args):
		env = list(self.environment)
		for element in args:
			if element not in env: return False
			env.remove(element)
		return len(env) == 0

def get_types(atoms, bonds):
	lst = [CustomAtom(symb) for symb in atoms.get_chemical_symbols()]
	
	for index, atom in enumerate(lst):
		atom.index = index
		for bond in bonds:
			n = None
			if bond[0].index == index:
				n = lst[bond[1].index]
			elif bond[1].index == index:
				n = lst[bond[0].index]
			if n:
				atom.add_bonded(n)
	
	for atom in lst:
		atom.rings = find_rings(atom)
	
	result = []
	for atom in lst:
		type = eval(atom.symbol).__call__(atom)
		if type == None:
			raise Exception('No CHARMM type for %s, %s' % (atom.symbol, atom.environment))
		else: 
			#print type[0], ':', type[1]
			result.append(type[0])
	return result
	
def find_rings(atom):
	rings = []
	def recursive_finder(chain):
		for n in chain[-1].bonded_atoms:
			if len(chain) > 2 and n == chain[0]: rings.append(chain)
			elif not n in chain: recursive_finder(chain + [n])
			
	recursive_finder([atom])
	return rings

# Charmm General force field types

# Utility function
def is_carbonyl_carbon(atom):
	if atom.symbol != 'C': return False
	for n in atom.bonded_atoms:
		if n.symbol == 'O'and len(n.environment) == 1:
			return True
	return False

## Carbon
def C(atom):
	neighs = atom.environment
	if len(neighs) == 4:
		if neighs.count('H') == 4:
			return 'CG331', 'aliphatic C for methyl group (-CH3)'
		elif neighs.count('H') == 3:
			return 'CG331', 'aliphatic C for methyl group (-CH3)'
		elif neighs.count('H') == 2:
			return 'CG321', 'aliphatic C for CH2'
		elif neighs.count('H') == 1:
			return 'CG311', 'aliphatic C with 1 H, CH'

	elif len(neighs) == 3:
		if atom.in_6_mem_ring():
			ring = atom.rings[0]
			nitrogen_neighs = [n for n in atom.bonded_atoms if n.symbol == 'N']
			if atom.in_5_mem_ring():
				return 'CG2RC0', '6/5-mem ring bridging C, guanine C4,C5, trp'
			elif len(nitrogen_neighs) > 1 and \
			   any([at for at in nitrogen_neighs if len(at.environment) == 2]):
				return 'CG2R64', '6-mem aromatic amidine and guanidine carbon (between 2 or 3 Ns and double-bound to one of them), NA, PYRM'
			elif is_carbonyl_carbon(atom):
				return 'CG2R63', '6-mem aromatic amide carbon (NA) (and other 6-mem aromatic carbonyls?)'
			elif any([is_carbonyl_carbon(member) for member in ring]):
				return 'CG2R62', '6-mem aromatic C for protonated pyridine (NIC) and rings containing carbonyls (see CG2R63) (NA)'
			else:
				return 'CG2R61', '6-mem aromatic C'
		elif atom.in_5_mem_ring():
			for n in atom.neighbors_of_type('N'):
				if len(n.bonded_atoms) == 2:
					if len([at for at in atom.bonded_atoms if at.heteroatom()]) > 1:
						return 'CG2R53', '5-mem ring, double bound to N and adjacent to another heteroatom, purine C8, his CE1 (0,+1), 2PDO, kevo'
					else:
						return 'CG2R52', '5-mem ring, double bound to N, PYRZ, pyrazole'
			# Else (not double bound to N)
			return 'CG2R51', '5-mem ring, his CG, CD2(0), trp'
		
		elif is_carbonyl_carbon(atom):
			if 'N' in neighs:
				return 'CG2O1', 'carbonyl C: amides'
			elif neighs.count('O') == 2:
				return 'CG2O2', 'carbonyl C: esters, [neutral] carboxylic acids'
			elif 'H' in neighs:
				return 'CG2O4', 'carbonyl C: aldehydes'
			else:
				return 'CG2O5', 'carbonyl C: ketones'
		
		elif neighs.count('H') == 1:
			return 'CG2D1', 'alkene; RHC= ; imine C'
		elif neighs.count('H') == 2:
			return 'CG2D2', 'alkene; H2C='
	
	elif len(neighs) == 2:
		return 'CG1T1', 'alkyn C'

# Hydrogen
def H(atom):
	if atom.environment[0] == 'C':
		carbon = atom.bonded_atoms[0]
		c_neighs = carbon.environment
		if carbon.in_6_mem_ring():
			if any([at for at in carbon.environment if not at in ('C', 'H')]) or \
			   any([at for at in carbon.rings[0] if is_carbonyl_carbon(at)]):
				return 'HGR62', 'nonnpolar H, neutral 6-mem planar ring C adjacent to heteroatom'
				       # or the ring contains carbonyl
			else:
				return 'HGR61', 'aromatic H'
		elif carbon.in_5_mem_ring():
			if any([at for at in carbon.environment if not at in ('C', 'H')]):
				return 'HGR52', 'Aldehyde H, formamide H (RCOH); nonpolar H, neutral 5-mem planar ring C adjacent to heteroatom or + charge'
			else:
				return 'HGR51', 'aromatic H'
		elif len(c_neighs) == 4:
			if 'N' in c_neighs:
				for n in carbon.neighbors_of_type('N'):
					if n.has_environment('C', 'H', 'H'):
						return 'HGAAM2', 'aliphatic H, NEUTRAL methylamine'
					elif n.has_environment('C', 'C', 'H'):
						return 'HGAAM1', 'aliphatic H, NEUTRAL dimethylamine'
					elif n.has_environment('C', 'C', 'C'):
						return 'HGAAM0', 'aliphatic H, NEUTRAL trimethylamine'
			
			elif c_neighs.count('H') == 1:
				return 'HGA1', 'alphatic proton, CH'
			elif c_neighs.count('H') == 2:
				return 'HGA2', 'alphatic proton, CH2'
			elif c_neighs.count('H') == 3:
				return 'HGA3', 'alphatic proton, CH3'
			elif c_neighs.count('H') == 4:
				return 'HGA3', 'alphatic proton, CH3'
		elif len(c_neighs) == 3:
			if is_carbonyl_carbon(carbon):
				return 'HGR52', 'Aldehyde H, formamide H (RCOH); nonpolar H, neutral 5-mem planar ring C adjacent to heteroatom or + charge'
			elif c_neighs.count('H') == 1:
				return 'HGA4', 'alkene proton; RHC='
			elif c_neighs.count('H') == 2:
				return 'HGA5', 'alkene proton; H2C=CR'
		elif len(c_neighs) == 2:
			# What to use here?
			return 'HGA1', 'alphatic proton, CH'
	
	elif atom.has_environment('O'):
		if list(atom.bonded_atoms[0].environment) == ['H', 'H']:
			return 'HGTIP3', 'polar H, TIPS3P WATER HYDROGEN'
		else:
			return 'HGP1', 'polar H'
	
	elif atom.has_environment('N'):
		nitrogen = atom.bonded_atoms[0]
		n_neighs = nitrogen.environment
		if nitrogen.aromatic():
			return 'HGP1', 'polar H'
		elif nitrogen.has_environment('C', 'H', 'H'):
			carbon = nitrogen.neighbors_of_type('C')[0]
			if carbon.n_neigbors() == 4:
				return 'HGPAM2', 'polar H, NEUTRAL methylamine'
			elif set(carbon.environment) == set(['N', 'C']):
				return 'HGP4', 'polar H, neutral conjugated -NH2 group (NA bases)'
			else:
				return 'HGP1', 'polar H'
		elif nitrogen.has_environment('C', 'C', 'H'):
			return 'HGPAM1', 'polar H, NEUTRAL dimethylamine'
		elif nitrogen.has_environment('H', 'H', 'H'):
			return 'HGPAM3', 'polar H, NEUTRAL ammonia'
		elif len(n_neighs) == 2:
			return 'HGP3', 'polar H, thiol'
		else:
			return 'HGP1', 'polar H'

# Oxygen
def O(atom):
	
	if atom.environment.count('H') == 1:
		return 'OG311', 'hydroxyl oxygen'
	
	if atom.environment.count('H') == 2:
		return 'OGTIP3', 'TIPS3P WATER OXYGEN'
	
	if atom.in_5_mem_ring():
		return 'OG2R50', 'FURA, furan'
	
	elif len(atom.environment) == 2:
		carbon_neighbors = [n for n in atom.bonded_atoms if n.symbol == 'C']
		for c in carbon_neighbors:
			# Is this an ester?
			oxygens = [n for n in carbon.bonded_atoms if n.symbol == 'O']
			if any([len(ox.environment) == 1] for ox in oxygens):
				return 'OG302', 'ester -O-'
			else:
				return 'OG301', 'ether -O-'
				
	
	elif list(atom.environment) == ['C']:
		carbon = atom.bonded_atoms[0]
		
		if list(carbon.environment) == ['O', 'O']:
			return 'OG2D5', 'CO2 oxygen'
		
		elif carbon.in_6_mem_ring():
			return 'OG2D4', '6-mem aromatic carbonyl oxygen (nucleic bases)'
		
		elif carbon.environment.count('C') == 2:
			return 'OG2D3', 'carbonyl O: ketones'
		else:
			return 'OG2D1', 'carbonyl O: amides, esters, [neutral] carboxylic acids, aldehydes, uera'

# Nitrogen
def N(atom):
	neighs = atom.environment 
	if list(neighs) == ['C']:
		return 'NG1T1', 'N for cyano group'

	elif atom.in_6_mem_ring():
		if len(neighs) == 2:
			for n in atom.bonded_atoms:
				if n.heteroatom() or \
				   len([at for at in n.bonded_atoms if at.heteroatom() and at.aromatic()]) > 1:
					return 'NG2R62', 'double bound 6-mem planar ring with heteroatoms in o or m, pyrd, pyrm'
			# Else
			return 'NG2R60', 'double bound neutral 6-mem planar ring, pyr1, pyzn'
		elif len(neighs) == 3:
			return 'NG2R61', 'single bound neutral 6-mem planar ring imino nitrogen; glycosyl linkage'
		
	elif atom.in_5_mem_ring():
		if len(neighs) == 2:
			return 'NG2R50', 'double bound neutral 5-mem planar ring, purine N7'
		elif len(neighs) == 3:
			return 'NG2R51', 'single bound neutral 5-mem planar (all atom types sp2) ring, his, trp pyrrole (fused)'
		
	elif len(neighs) == 3:
		if any([is_carbonyl_carbon(n) for n in atom.bonded_atoms]):
			if neighs.count('H') == 1:
				return 'NG2S1', 'peptide nitrogen (CO=NHR)'
			elif  neighs.count('H') == 2:
				return 'NG2S2', 'terminal amide nitrogen (CO=NH2)'
		
		elif neighs.count('C') == 3:
			return 'NG301', 'neutral trimethylamine nitrogen'
		elif neighs.count('C') == 2 and neighs.count('H') == 1:
			return 'NG311', 'neutral dimethylamine nitrogen'
		elif neighs.count('C') == 1 and neighs.count('H') == 2:
			for n in atom.bonded_atoms:
				if n.symbol == 'C' and n.aromatic():
					return 'NG2S3', 'external amine ring nitrogen (planar/aniline), phosphoramidate'
			# Else
			return 'NG321', 'neutral methylamine nitrogen'

		elif list(neighs) == ['H']*3:
			return 'NG331', 'neutral ammonia nitrogen'
	
	elif len(neighs) == 2:
		if 'C' in neighs:
			return 'NG2D1', "N for neutral imine/Schiff's base (C=N-R, acyclic amidine, gunaidine)"
	
from itertools import combinations
from ase.data import covalent_radii
def detectBonds(atoms):
	bonds = []
	tolerance = 1.2
	for a1, a2 in combinations(atoms, 2):
		bondLength = covalent_radii[a1.number] + covalent_radii[a2.number]
		if atoms.get_distance(a1.index, a2.index, mic=True) < bondLength*tolerance:
			bonds.append([a1.index, a2.index])
	return bonds

if __name__ == '__main__':
	from  ase.data import s22
	for name in s22.s22:
		print '\n***', name, '***'
		atoms = s22.create_s22_system(name)
		get_types(atoms, detectBonds(atoms))