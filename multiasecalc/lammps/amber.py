from lammpsbase import LAMMPSBase
import read_gromacs_params
import typing, ambertypes
import warnings
import os

class AMBER(LAMMPSBase):
	
	def __init__(self, gromacs_dir, label='amber', inner_cutoff = 8.0, outer_cutoff=10.0, **kwargs):
		LAMMPSBase.__init__(self, label, **kwargs)
		
		self.parameters.units          = 'real'
		self.parameters.pair_style     = 'lj/charmm/coul/charmm %f %f' % (inner_cutoff, outer_cutoff)
		self.parameters.bond_style     = 'harmonic'
		self.parameters.angle_style    = 'harmonic'
		self.parameters.dihedral_style = 'charmm'
		#self.parameters.improper_style = 'harmonic'
		self.parameters.special_bonds  = 'amber'
		
		# Read force field file
		ffnonbonded = os.path.join(gromacs_dir, 'amber03.ff', 'ffnonbonded.itp')
		ffbonded = os.path.join(gromacs_dir, 'amber03.ff', 'ffbonded.itp')
		
		self.ff_data = read_gromacs_params.read((ffnonbonded, ffbonded))
		self.type_resolver = typing.TypeResolver(ambertypes.data)
	
	def atom_types(self, atoms):
		self.type_resolver.resolve_atoms(atoms)
		return atoms.info['atom_types']
		
	def prepare_calculation(self, atoms, data):
		if (atoms.get_charges() == 0).all():
			warnings.warn("No partial charges set! There won't be any Coulomb interactions.")
