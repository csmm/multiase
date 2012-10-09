from lammpsbase import LAMMPSBase
import read_compass_params
import typing, compasstypes

class COMPASS(LAMMPSBase):
	
	def __init__(self, label='compass', ff_file_path='compass.frc', pair_cutoff=10.0, **kwargs):
		LAMMPSBase.__init__(self, label, **kwargs)
		
		self.parameters.units          = 'real'
		self.parameters.pair_style     = 'lj/class2/coul/cut %f' % pair_cutoff
		self.parameters.bond_style     = 'class2'
		self.parameters.angle_style    = 'class2'
		self.parameters.dihedral_style = 'class2'
		self.parameters.improper_style = 'class2'
		self.parameters.special_bonds  =  'lj/coul 0.0 0.0 1.0 dihedral yes'
		
		self.ff_data, self.bond_increments = read_compass_params.read(open(ff_file_path))
		self.type_resolver = typing.TypeResolver(compasstypes.data)
	
	def atom_types(self, atoms):
		try:
			return [self.type_resolver.resolve(atom).type for atom in atoms]
		except AttributeError:
			raise RuntimeError('Could not resolve all types!')
		
	
	def prepare_calculation(self, atoms, data):
		bonds = atoms.get_array('bonds')
		types = data.atom_types
		
		for i in range(len(atoms)):
			q = 0
			t1 = types[i]
			for j in bonds[i]:
				t2 = types[j]
				increment = self.bond_increments.get((t1,t2))
				if increment:
					q += increment[0]
				else:
					increment = self.bond_increments.get((t2,t1))
					if not increment:
						raise RuntimeError('No bond increment for (%s, %s)' % (t1, t2))
					q += increment[1]
			atoms[i].charge = q
		