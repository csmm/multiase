from lammpsbase import LAMMPSBase, FFData
from tempfile import NamedTemporaryFile
import numpy as np
import shutil, os

def get_element_order(ff_file):
	f = open(ff_file)
	
	order = []
	for line in f:
		fields = line.strip().split()
		if len(fields) == 9 and fields[0][0].isupper():
			order.append(fields[0])
	return order

class ReaxFF(LAMMPSBase):
	
	def __init__(self, ff_file_path, label='reaxff', specorder=None,
		implementation='C', update_charges=True, save_bond_orders=False, debug_energy=False, **kwargs):
		
		LAMMPSBase.__init__(self, label, update_charges = update_charges, **kwargs)
		
		
		self.parameters.atom_style = 'charge'
		self.parameters.pair_style = 'reax'
		if implementation == 'C':
			self.parameters.pair_style += '/c NULL'
			self.parameters.extra_cmds.append('fix qeq all qeq/reax 1 0.0 10.0 1.0e-6 reax/c')
		
		self.parameters.units      = 'real'
		
		if not specorder:
			self.specorder = get_element_order(ff_file_path)
		else:
			self.specorder = specorder
		
		self.ff_file = ff_file_path
		self.ff_data = FFData()
		
		self.save_bond_orders = save_bond_orders
		"""
		if save_bond_orders:
			self.bondinfo_file = NamedTemporaryFile(mode='r', prefix='bonds_'+label,
									dir=self.tmp_dir, delete=(not self.keep_tmp_files))
			self.parameters.extra_cmds.append('fix bondinfo all reax/bonds 1 %s\n' % (self.bondinfo_file.name))
		
		if debug_energy:
			energy_terms = 'eb', 'ea', 'elp', 'emol', 'ev', 'epen', 'ecoa', 'ehb', 'et', 'eco', 'ew', 'ep', 'efi', 'eqeq'
			self.parameters.extra_cmds += ['compute reax all pair reax']
			self.parameters.extra_cmds += ['variable %s equal c_reax[%i]' % (et, i+1) for i, et in enumerate(energy_terms)]
			
			self._custom_thermo_args += ['v_%s' % et for et in energy_terms]
		"""
	
	def atom_types(self, atoms):
		return atoms.get_chemical_symbols()
	
	def prepare_calculation(self, atoms, data):
		typeorder = data.atom_typeorder
		element_indices = [str(self.specorder.index(el)+1) for el in typeorder]
		self.parameters.pair_coeffs = ['* * %s ' % (self.ff_file) + ' '.join(element_indices)]
	
	
	def set_dumpfreq(self, freq):
		LAMMPSBase.set_dumpfreq(self, freq)
		#if self.save_bond_orders:
		#	self.lammps_input.write('fix bondinfo all reax/bonds %s %s\n' % (freq, self.bondinfo_file.name))
	
	def read_lammps_trj(self, *args, **kwargs):
		LAMMPSBase.read_lammps_trj(self, *args, **kwargs)
		#if self.save_bond_orders:
		#	self.read_bondinfo()
		
	def read_bondinfo(self):
		f = self.bondinfo_file
		
		bondarr = [[]]*len(self.atoms)
		orderarr = [[]]*len(self.atoms)
		aboarr = np.empty([len(self.atoms)])
		nlparr = np.empty([len(self.atoms)])
		
		for line in f:
			if '#' in line: continue
			fields = line.strip().split()
			index = int(fields[0])-1
			nbonds = int(fields[2])
			bonds = [int(field)-1 for field in fields[3:3+nbonds]]
			bondarr[index] = bonds
			startind = 3+nbonds+1
			orders = [float(field) for field in fields[startind:startind+nbonds]]
			orderarr[index] = orders
			
			aboarr[index] = fields[startind+nbonds]
			nlparr[index] = fields[startind+nbonds+1]
			
		bondarr = np.array(bondarr, dtype=object)
		orderarr = np.array(orderarr, dtype=object)
		
		atoms = self.atoms
		atoms.set_array('bonds', bondarr)
		atoms.set_array('bond_orders', orderarr)
		atoms.set_array('abos', aboarr)
		atoms.set_array('nlps', nlparr)