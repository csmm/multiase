from lammpsbase import LAMMPSBase
from ffdata import SequenceType
import read_gromacs_params
import typing, oplsaatypes
import warnings
import os

class OPLSAA(LAMMPSBase):
	
	def __init__(self, gromacs_dir, label='opls-aa', pair_cutoff=10.0, kspace=False, **kwargs):
		LAMMPSBase.__init__(self, label, **kwargs)
		
		self.parameters.units          = 'real'
		if kspace:
			self.parameters.pair_style     = 'lj/cut/coul/long %f' % pair_cutoff
			self.parameters.kspace_style   = 'ewald 0.0001'
		else:
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
		
		# Add parameters for flexible TIP3P water
		self.ff_data.add('bond', SequenceType(['OW', 'HW']), 'Bond Coeffs', [450, 0.9572])
		self.ff_data.add('angle', SequenceType(['HW', 'OW', 'HW']), 'Angle Coeffs', [55, 104.52])
		
	def atom_types(self, atoms):
		self.type_resolver.resolve_atoms(atoms)
		return atoms.info['atom_types']
	
	def set_charges(self, atoms, atom_types):
		charges = atoms.get_charges()
		for i in range(len(atoms)):
			charges[i] = self.charges[atom_types[i]]
		atoms.set_charges(charges)