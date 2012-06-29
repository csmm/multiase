from lammpsBase import LAMMPSBase
from ase.data import covalent_radii
from ase.atoms import Atoms
from ase.atom import Atom, atomproperty
import ase.atom
from itertools import combinations, product, permutations
import readFrc, compassTypes
import numpy as np

class COMPASS(LAMMPSBase):
	
	def __init__(self, label='compass', ff_file_path='compass.frc', pair_cutoff=10.0, **kwargs):
		LAMMPSBase.__init__(self, label, **kwargs)
		
		self.parameters.units          = 'real'
		self.parameters.pair_style     = 'lj/class2/coul/cut %f' % pair_cutoff
		self.parameters.bond_style     = 'class2'
		self.parameters.angle_style    = 'class2'
		self.parameters.dihedral_style = 'class2'
		self.parameters.improper_style = 'class2'
		
		# Read force field file
		self.use_ff_file(ff_file_path)
		self.ff_parameters = readFrc.read(open(ff_file_path))
		
		# TODO
		# pair_modify     tail yes
		# special_bonds   lj/coul 0.0 0.0 1.0 dihedral yes	
	
	def set_bonds(self, bonds):
		self.raw_bonds = [(bond[0].index, bond[1].index) for bond in bonds]
	
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
		
		# Add a compass type property to ASE atoms
		if not atoms.has('compass_types'):
			atoms.new_array('compass_types', 
				np.array([compassTypes.getType(a, self.atomNeighbors(a)) for a in atoms]))
		
		ase.atom.names['compass_type'] = ('compass_types', '')
		ase.atom.Atom.compass_type = ase.atom.atomproperty('compass_type', 'COMPASS atom type')
		
		
		# Calculate charges
		for atom in atoms:
			atom.charge = self.calculateCharge(atom,
				ffParameters.bondIncrement, ffParameters.equivalence)
		
		# Set types for bonds, angles, etc
		for bond in bonds:
			bond.compass_type = self.standardOrder(tuple(atom.compass_type for atom in bond))
			
		for angle in angles:
			angle.compass_type = self.standardOrder(tuple(atom.compass_type for atom in angle))
			
		for dihedral in dihedrals:
			dihedral.compass_type = self.standardOrder(tuple(atom.compass_type for atom in dihedral))
		
		# An improper ijkl consist of bonds (j,i) (j,k) (j,l)
		for improper in impropers:
			types = [atom.compass_type for atom in improper]
			i, k, l = sorted((types[0], types[2], types[3]))
			improper.compass_type = (i, types[1], k, l)
		
		# Make the sets ordered
		atomTypes     = list(set(a.compass_type for a in atoms))
		bondTypes     = list(set(b.compass_type for b in bonds))
		angleTypes    = list(set(a.compass_type for a in angles))
		dihedralTypes = list(set(d.compass_type for d in dihedrals))
		improperTypes = list(set(i.compass_type for i in impropers))
	
		def build_coeff_table(title, formalTypes, actualInteraction, dataTable):
			data = []
			for index, type in enumerate(formalTypes, 1):
				actualType = tuple(ffParameters.equivalence[ft][actualInteraction] for ft in type)
				if not actualType in dataTable:
					actualType = actualType[::-1]
				coeffs = dataTable.get(actualType)
				if coeffs == None:
					print 'Warning: no {0} coeffs for {1}'.format(title, actualType)
					coeffs = dataTable['empty']
				data.append(coeffs)
			return dict(title=title, data=data)
			
		argument_sets = [
			# title                formalTypes  actualInteraction  dataTable
			('Pair Coeffs', [(at,) for at in atomTypes], 'nonbonded', ffParameters.pair),
			('Bond Coeffs',           bondTypes, 'bonded',      ffParameters.bond),
			('Angle Coeffs',          angleTypes, 'angle',      ffParameters.angle),
			('BondAngle Coeffs',      angleTypes, 'bonded',     ffParameters.bondAngle),
			('BondBond Coeffs',       angleTypes, 'bonded',     ffParameters.bondBond),
			('Dihedral Coeffs',       dihedralTypes, 'torsion', ffParameters.torsion),
			('EndBondTorsion Coeffs', dihedralTypes, 'bonded',  ffParameters.endBondTorsion),
			('MiddleBondTorsion Coeffs', dihedralTypes, 'bonded', ffParameters.middleBondTorsion),
			('BondBond13 Coeffs',     dihedralTypes, 'bonded',  ffParameters.bondBond13),
			('AngleTorsion Coeffs',   dihedralTypes, 'angle',   ffParameters.angleTorsion),
			('AngleAngleTorsion Coeffs', dihedralTypes, 'angle', ffParameters.angleAngleTorsion)
			]
		
		# Populate self.data.coeff_tables
		for argument_set in argument_sets:
			self.data.coeff_tables.append(build_coeff_table(*argument_set))
			
			
		# Impropers have a different symmetry
		def build_improper_coeffs(title, formalTypes, actualInteraction, dataTable):
			data = []
			for type in improperTypes:
				actualType = tuple(ffParameters.equivalence[ft][actualInteraction] for ft in type)
				coeffs = None
				for i, k, l in permutations((actualType[0], actualType[2], actualType[3])):
					coeffs = dataTable.get((i, actualType[1], k, l))
					if coeffs:
						break
				if not coeffs:
					print 'Warning: no {0} coeffs for {1}'.format(title, actualType)
					coeffs = dataTable['empty']
				data.append(coeffs)
			return dict(title=title, data=data)
				
		argument_sets = [
			('Improper Coeffs', improperTypes, 'torsion', ffParameters.improper),
			('AngleAngle Coeffs', improperTypes, 'angle', ffParameters.angleAngle)
			]
		
		# Populate self.data.coeff_tables
		for argument_set in argument_sets:
			self.data.coeff_tables.append(build_improper_coeffs(*argument_set))
			
		self.data.atom_types = [
			atomTypes.index(a.compass_type)+1 for a in atoms]
			
		self.data.bonds = []
		for bond in bonds:
			type = bondTypes.index(bond.compass_type)+1
			self.data.bonds.append([type] + [atom.index+1 for atom in bond])
			
		self.data.angles = []
		for angle in angles:
			type = angleTypes.index(angle.compass_type)+1
			self.data.angles.append([type] + [atom.index+1 for atom in angle])
		
		self.data.dihedrals = []
		for dihedral in dihedrals:
			type = dihedralTypes.index(dihedral.compass_type)+1
			self.data.dihedrals.append([type] + [atom.index+1 for atom in dihedral])
			
		self.data.impropers = []
		for improper in impropers:
			type = improperTypes.index(improper.compass_type)+1
			self.data.impropers.append([type] + [atom.index+1 for atom in improper])
		
		# Masses
		massDict = dict((a.compass_type, a.mass) for a in atoms)
		self.data.masses = [massDict[tp] for tp in atomTypes]
		
		# Clean up
		del Atom.compass_type
		del ase.atom.names['compass_type']
		
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
	"""
	def detectBonds(self):
		bonds = []
		tolerance = 1.2
		for a1, a2 in combinations(self.atoms, 2):
			bondLength = covalent_radii[a1.number] + covalent_radii[a2.number]
			if self.atoms.get_distance(a1.index, a2.index, mic=True) < bondLength*tolerance:
				bonds.append(Bond((a1, a2)))
		return bonds
	"""
	
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
				impr = Improper((triplet[0], atom, triplet[1], triplet[2]))
				impropers.append(impr)
		return impropers
	
	
	def calculateCharge(self, atom, bondIncrements, equivalence):
		q = 0
		neighs = self.atomNeighbors(atom)
		atomType = equivalence[atom.compass_type]['bonded']
		neighTypes = [equivalence[n.compass_type]['bonded'] for n in neighs]
		for pair in ((atomType, nt) for nt in neighTypes):
			if pair in bondIncrements:
				q += bondIncrements[pair][0]
			elif pair[::-1] in bondIncrements:
				q -= bondIncrements[pair[::-1]][0]
			else:
				print 'Warning: no bond increment for ', pair
		return q
	
	
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
