from lammpsbase import LAMMPSBase
from ase.data import covalent_radii
from ase.atoms import Atoms
from ase.atom import Atom, atomproperty
import ase.atom
from itertools import combinations, product, permutations
import read_charmm_params, charmm_types
import numpy as np

class CHARMM(LAMMPSBase):
	
	def __init__(self, label='charmm', ff_file_path=None, pair_cutoff=10.0, **kwargs):
		LAMMPSBase.__init__(self, label, **kwargs)
		
		self.parameters.units          = 'real'
		self.parameters.pair_style     = 'lj/charmm/coul/charmm 8.0 10.0'
		self.parameters.bond_style     = 'harmonic'
		self.parameters.angle_style    = 'charmm'
		self.parameters.dihedral_style = 'charmm'
		self.parameters.improper_style = 'harmonic'
		self.parameters.special_bonds  = 'charmm'
		
		# Read force field file
		self.use_ff_file(ff_file_path)
		self.ff_parameters = read_charmm_params.read(open(ff_file_path))
	
	def prepare_calculation(self):
		""" Fill self.data """
		
		self.data.clear()
		
		atoms = self.atoms
		
		if not atoms.has('bonds'):
			atoms.set_array('bonds', self.detectBonds())
			
		bondarr = atoms.get_array('bonds')
		bonds = [Bond((a, b)) for a,b in combinations(atoms, 2) if b.index in bondarr[a.index]]
		self.bonds = bonds
		angles    = self.findAngles()
		dihedrals = self.findDihedrals()
		impropers = self.findImpropers()
		
		ffParameters = self.ff_parameters
		
		# Add a charmm type property to ASE atoms
		charmm_atom_types = np.array(charmm_types.get_types(atoms, bonds))
		if not atoms.has('charmm_types'):
			atoms.set_array('charmm_types', charmm_atom_types)
		
		ase.atom.names['charmm_type'] = ('charmm_types', '')
		ase.atom.Atom.charmm_type = ase.atom.atomproperty('charmm_type', 'CHARMM atom type')
		
		# Set types for bonds, angles, etc
		for bond in bonds:
			bond.charmm_type = self.standardOrder(tuple(atom.charmm_type for atom in bond))
			
		for angle in angles:
			angle.charmm_type = self.standardOrder(tuple(atom.charmm_type for atom in angle))
			
		for dihedral in dihedrals:
			dihedral.charmm_type = self.standardOrder(tuple(atom.charmm_type for atom in dihedral))
		
		# An improper ijkl consist of bonds (i,j) (i,k) (i,l)
		for improper in impropers:
			types = [atom.charmm_type for atom in improper]
			improper.charmm_type = (types[0],) + tuple(sorted(types[1:]))
		
		# Make the sets ordered
		atom_types     = list(set(a.charmm_type for a in atoms))
		bond_types     = list(set(b.charmm_type for b in bonds))
		angle_types    = list(set(a.charmm_type for a in angles))
		dihedral_types = list(set(d.charmm_type for d in dihedrals))
		improper_types = list(set(i.charmm_type for i in impropers))
		
		def build_coeff_table(title, types, dataTable):
			data = []
			for index, type in enumerate(types, 1):
				if not type in dataTable:
					type = type[::-1]
				coeffs = dataTable.get(type)
				if coeffs == None:
					print 'Warning: no {0} for {1}'.format(title, type)
					coeffs = dataTable['empty']
				data.append(coeffs)
			return dict(title=title, data=data)
		
		argument_sets = [
			# title                formalTypes  actualInteraction  dataTable
			('Pair Coeffs', [(at,) for at in atom_types], ffParameters.nonbonded),
			('Bond Coeffs',           bond_types,  ffParameters.bond),
			('Angle Coeffs',          angle_types,  ffParameters.angle),
			('Dihedral Coeffs',       dihedral_types, ffParameters.dihedral),
			]
		
		# Populate self.data.coeff_tables
		for argument_set in argument_sets:
			self.data.coeff_tables.append(build_coeff_table(*argument_set))
	
		# Build improper table
		data = []
		for type in improper_types:
			coeffs = None
			for j, k, l in permutations(type[1:]):
				coeffs = ffParameters.improper.get((type[0], j, k, l))
				if coeffs:
					#print 'Found improper coeffs for %s' % (type,)
					break
			if not coeffs:
				#print 'Warning: no improper coeffs for %s' % (type,)
				coeffs = ffParameters.improper['empty']
			data.append(coeffs)
		self.data.coeff_tables.append(dict(title='Improper Coeffs', data=data))
		
		
		self.data.atom_types = [
			atom_types.index(a.charmm_type)+1 for a in atoms]
			
		self.data.bonds = []
		for bond in bonds:
			type = bond_types.index(bond.charmm_type)+1
			self.data.bonds.append([type] + [atom.index+1 for atom in bond])
			
		self.data.angles = []
		for angle in angles:
			type = angle_types.index(angle.charmm_type)+1
			self.data.angles.append([type] + [atom.index+1 for atom in angle])
		
		self.data.dihedrals = []
		for dihedral in dihedrals:
			type = dihedral_types.index(dihedral.charmm_type)+1
			self.data.dihedrals.append([type] + [atom.index+1 for atom in dihedral])
			
		self.data.impropers = []
		for improper in impropers:
			type = improper_types.index(improper.charmm_type)+1
			self.data.impropers.append([type] + [atom.index+1 for atom in improper])
		
		# Masses
		massDict = dict((a.charmm_type, a.mass) for a in atoms)
		self.data.masses = [massDict[tp] for tp in atom_types]
			
		# Clean up
		del Atom.charmm_type
		del ase.atom.names['charmm_type']
	
	def detectBonds(self):
		from ase.data import covalent_radii
		atoms = self.atoms
		tolerance = 1.15
		def bonded(i, j):
			bondLength = covalent_radii[atoms[i].number] + covalent_radii[atoms[j].number]
			return atoms.get_distance(i, j, mic=True) < bondLength*tolerance
				
		N = len(atoms)
		bonds = np.empty((N), dtype=object)
		for i in range(N):
			bonds[i] = [j for j in range(N) if bonded(i, j)]
		return bonds
		
	
	def atomNeighbors(self, atom):
		return [b.otherAtom(atom) for b in self.bonds if atom in b]
	
	
	def findAngles(self):
		angles = []
		for atom in self.atoms:
			neighs = self.atomNeighbors(atom)
			for pair in combinations(neighs, 2):
				angles.append(Angle((pair[0], atom, pair[1])))
		return angles
				
	def findDihedrals(self):
		dihedrals = []
		for bond in self.bonds:
			neighbors1 = filter(lambda b: bond[0] in b and b != bond, self.bonds)
			neighbors2 = filter(lambda b: bond[1] in b and b != bond, self.bonds)
			for b1, b2 in product(neighbors1, neighbors2):
				first = b1.otherAtom(bond[0])
				last  = b2.otherAtom(bond[1])
				dh = Dihedral((first, bond[0], bond[1], last))
				dihedrals.append(dh)
		return dihedrals

		
	def findImpropers(self):
		impropers = []
		for atom in self.atoms:
			neighs = self.atomNeighbors(atom)
			for triplet in combinations(neighs, 3):
				# In CHARMM improper consists of bonds ij, ik, il
				impr = Improper((atom, triplet[0], triplet[1], triplet[2]))
				impropers.append(impr)
		return impropers
	
	
	
	def standardOrder(self, sequence):
		""" Order a sequence of atom types in a predictable way """
		if sequence[0] > sequence[-1]: sequence = sequence[::-1]
		return sequence
		
			
class Bond(tuple):
	def otherAtom(self, atom):
		if   self[0].index == atom.index: return self[1]
		elif self[1].index == atom.index: return self[0]
		
	def __contains__(self, atom):
		return atom.index == self[0].index or atom.index == self[1].index
		
class Angle(tuple):
	pass

class Dihedral(tuple):
	pass

class Improper(tuple):
	pass
 
