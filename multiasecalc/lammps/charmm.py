from lammpsbase import LAMMPSBase
import read_charmm_params
import typing, charmmtypes

class CHARMM(LAMMPSBase):
	
	def __init__(self, ff_file_path, label='charmm', pair_cutoff=10.0, **kwargs):
		LAMMPSBase.__init__(self, label, **kwargs)
		
		self.parameters.units          = 'real'
		self.parameters.pair_style     = 'lj/charmm/coul/charmm 8.0 10.0'
		self.parameters.bond_style     = 'harmonic'
		self.parameters.angle_style    = 'charmm'
		self.parameters.dihedral_style = 'charmm'
		self.parameters.improper_style = 'harmonic'
		self.parameters.special_bonds  = 'charmm'
		
		# Read force field file
		self.ff_data = read_charmm_params.read(open(ff_file_path))
		self.type_resolver = typing.TypeResolver(charmmtypes.data)
	
	def atom_types(self, atoms):
		return [self.type_resolver.resolve(atom).type for atom in atoms]
		