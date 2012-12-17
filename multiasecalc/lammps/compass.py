from lammpsbase import LAMMPSBase
from ffdata import SequenceType
import read_compass_params
import typing, compasstypes

class COMPASS(LAMMPSBase):
	
	def __init__(self, ff_file_path, label='compass', pair_cutoff=9.5, kspace=False, **kwargs):
		LAMMPSBase.__init__(self, label, **kwargs)
		
		self.parameters.units          = 'real'
		
		if kspace:
			self.parameters.pair_style     = 'lj/class2/coul/long %f' % pair_cutoff
			self.parameters.kspace_style   = 'ewald 0.0001'
		else:
			self.parameters.pair_style     = 'lj/class2/coul/cut %f' % pair_cutoff
			
		self.parameters.bond_style     = 'class2'
		self.parameters.angle_style    = 'class2'
		self.parameters.dihedral_style = 'class2'
		self.parameters.improper_style = 'class2'
		self.parameters.special_bonds  =  'lj/coul 0.0 0.0 1.0 dihedral yes'
		
		self.ff_data, self.bond_increments = read_compass_params.read(open(ff_file_path))
		self.type_resolver = typing.TypeResolver(compasstypes.data)
	
	def atom_types(self, atoms):
		self.type_resolver.resolve_atoms(atoms)
		return atoms.info['atom_types']
	
	def set_charges(self, atoms, atom_types):
		bonds = atoms.info['bonds']
		eq = self.ff_data.equivalence
		
		for i in range(len(atoms)):
			q = 0
			t1 = eq[atom_types[i]]['bond']
			for j in bonds[i]:
				t2 = eq[atom_types[j]]['bond']
				increment = self.bond_increments.get((t1, t2))
				if increment:
					q += increment[0]
				else:
					increment = self.bond_increments.get((t2, t1))
					if not increment:
						raise RuntimeError('No bond increment for (%s, %s)' % (t1, t2))
					q += increment[1]
			atoms[i].charge = q
		