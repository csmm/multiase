from lammpsbase import LAMMPSBase
import read_gromacs_params
import typing, oplsaatypes
import warnings
import os

class OPLSAA(LAMMPSBase):
	
	def __init__(self, gromacs_dir, label='opls-aa', pair_cutoff=10.0, **kwargs):
		LAMMPSBase.__init__(self, label, **kwargs)
		
		self.parameters.units          = 'real'
		self.parameters.pair_style     = 'lj/cut/coul/cut %f' % pair_cutoff
		self.parameters.pair_modify    = 'mix geometric'
		self.parameters.bond_style     = 'harmonic'
		self.parameters.angle_style    = 'harmonic'
		self.parameters.dihedral_style = 'multi/harmonic'
		#self.parameters.improper_style = 'harmonic'
		self.parameters.special_bonds  = 'lj/coul 0.0 0.0 0.5'
		
		# Read force field file
		ffnonbonded = os.path.join(gromacs_dir, 'oplsaa.ff', 'ffnonbonded.itp')
		ffbonded = os.path.join(gromacs_dir, 'oplsaa.ff', 'ffbonded.itp')
		
		self.ff_data, self.charges = read_gromacs_params.read((ffnonbonded, ffbonded), return_charges=True)
		self.type_resolver = typing.TypeResolver(oplsaatypes.data)
		
	def atom_types(self, atoms):
		self.type_resolver.resolve_atoms(atoms)
		return atoms.info['atom_types']
	
	def set_charges(self, atoms, atom_types):
		charges = atoms.get_charges()
		for i in range(len(atoms)):
			charges[i] = self.charges[atom_types[i]]
		atoms.set_charges(charges)