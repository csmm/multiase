from ase.calculators.interface import Calculator
import os, sys
import shutil, shlex
from subprocess import Popen, PIPE
from tempfile import mkdtemp, NamedTemporaryFile
import numpy as np
import unitconversion
from bonds import Bonds
from ffdata import FFData, SequenceType, ImproperType
from multiasecalc.lammps.io.lammps import read_lammps_dump

try:
	from itertools import combinations, permutations
except:
	# itertools.combinations and itertools.permutations are not in python 2.4
	def permutations(iterable, r=None):
		# permutations('ABCD', 2) --> AB AC AD BA BC BD CA CB CD DA DB DC
		# permutations(range(3)) --> 012 021 102 120 201 210
		pool = tuple(iterable)
		n = len(pool)
		if r is None: r = n
		if r > n:
			return
		indices = range(n)
		cycles = range(n, n-r, -1)
		yield tuple(pool[i] for i in indices[:r])
		while n:
			for i in reversed(range(r)):
				cycles[i] -= 1
				if cycles[i] == 0:
					indices[i:] = indices[i+1:] + indices[i:i+1]
					cycles[i] = n - i
				else:
					j = cycles[i]
					indices[i], indices[-j] = indices[-j], indices[i]
					yield tuple(pool[i] for i in indices[:r])
					break
			else:
				return
				
	def combinations(iterable, r):
		pool = tuple(iterable)
		n = len(pool)
		for indices in permutations(range(n), r):
			if sorted(indices) == list(indices):
				yield tuple(pool[i] for i in indices)


__ALL__ = ['LAMMPSBase']

# "End mark" used to indicate that the calculation is done
CALCULATION_END_MARK = '__end_of_ase_invoked_calculation__'

class LAMMPSParameters:
	def __init__(self, **pars):
		self.atom_style     = pars.get('atom_style', 'full')
		self.units          = pars.get('units')
		self.neighbor       = pars.get('neighbor')
		self.newton         = pars.get('newton')
		
		self.pair_style     = pars.get('pair_style')
		self.bond_style     = pars.get('bond_style')
		self.angle_style    = pars.get('angle_style')        
		self.dihedral_style = pars.get('dihedral_style')
		self.improper_style = pars.get('improper_style')
		self.kspace_style   = pars.get('kspace_style')
		self.special_bonds  = pars.get('special_bonds')
		self.pair_modify    = pars.get('pair_modify')
		
		# Pair coeffs to be specified in the input file
		self.pair_coeffs    = pars.get('pair_coeffs', [])
		
		self.extra_cmds     = pars.get('extra_cmds', [])
		
		# Atoms groups: tuples (group_id, list_of_indices)
		self.groups         = pars.get('groups', [])

class LAMMPSData:
	def __init__(self):
		self.clear()
		
	def clear(self):
		self.tables = None
		
		self.atom_types = None
		self.bonds = None
		self.angles = None
		self.dihedrals = None
		self.impropers = None
		
		self.atom_typeorder = None
		self.bond_typeorder = None
		self.angle_typeorder = None
		self.dihedral_typeorder = None
		self.improper_typeorder = None


class LAMMPSBase(Calculator):

	def __init__(self, label='lammps', tmp_dir=None, parameters={}, 
				update_charges = False, lammps_command=None,
				keep_alive=False, debug=False):
		"""The LAMMPS calculators object """

		self.label = label
		self.parameters = LAMMPSParameters(**parameters)
		self.data = LAMMPSData()
		self.ff_data = None
		self.forces = None
		self.atoms = None
		self.atoms_after_last_calc = None
		self.update_charges = update_charges
		self.keep_alive = keep_alive
		self.lammps_command = lammps_command
		self.debug = debug           
		self.lammps_process = LammpsProcess(log=debug)
		self.calls = 0
		
		self._custom_thermo_args = [
			'step', 'temp', 'press', 'cpu', 
			'pxx', 'pyy', 'pzz', 'pxy', 'pxz', 'pyz',
			'ke', 'pe', 'etotal',
			'vol', 'lx', 'ly', 'lz', 'atoms']
		
		self._dump_fields = [
			'id', 'type', 'x', 'y', 'z', 'vx',
			'vy', 'vz', 'fx', 'fy', 'fz', 'q']
		
		if tmp_dir is None:
			self.tmp_dir = mkdtemp(prefix='LAMMPS-')
		else:
			# If tmp_dir is pointing somewhere, don't remove stuff!
			self.debug = True
			self.tmp_dir=os.path.realpath(tmp_dir)
			if not os.path.isdir(self.tmp_dir):
				os.mkdir(self.tmp_dir, 0755)
		
		if debug:
			print 'LAMMPS (label: %s) running at %s' % (self.label, self.tmp_dir)
	
	def __del__(self):
		if not self.debug:
			shutil.rmtree(self.tmp_dir)
	
	def atom_types(self, atoms):
		""" Implement this method in subclasses"""
		raise NotImplementedException()
	
	def set_charges(self, atoms, atom_types):
		""" Implement this method in subclasses if needed"""
		pass
	
	def prepare_calculation(self, atoms, data):
		""" Implement this method in subclasses if needed"""
		pass

	def get_potential_energy(self, atoms):
		self.update(atoms)
		energy = self.lammps_process.get_thermo('pe')
		return self.to_ase_units(energy, 'energy')

	def get_forces(self, atoms):
		self.update(atoms)
		return self.forces

	def get_stress(self, atoms):
		self.update(atoms)
		thermo = self.lammps_process.get_thermo
		stress = [thermo(key) for key in ('pxx','pyy','pzz', 'pyz','pxz','pxy')]
		return self.to_ase_units(np.array(stress), 'stress')  
	
	def calculation_required(self, atoms, quantities=None):
		return atoms != self.atoms_after_last_calc
		
	def update(self, atoms):
		if not self.calculation_required(atoms): return
		self.setup_calculation(atoms)
		self.evaluate_forces()
		self.close_calculation()
	
	def minimize(self, atoms, etol=0, ftol=0, maxeval=100000, maxiter=100000000, min_style=None, relax_cell=False):
		if etol == 0 and ftol == 0: 
			raise RuntimeError('Specify at least one tolerance value!')
		ftol = self.from_ase_units(ftol, 'force')
		minimize_params = '%g %g %i %i' % (etol, ftol, maxeval, maxiter)
		
		self.setup_calculation(atoms)
		self.evaluate_forces()
		
		self.set_dumpfreq(999999)
		f = self.lammps_process
		if relax_cell:
			f.write('fix relax_cell all box/relax tri 0.0 nreset 20\n')
		if min_style:
			f.write('min_style %s\n' % min_style)
		f.write('minimize %s\n'  % minimize_params)
		if relax_cell:
			f.write('unfix relax_cell\n')
		f.write('print "%s"\n' % CALCULATION_END_MARK)
		f.flush()
		self.lammps_process.read_lammps_output()
		self.read_lammps_trj(update_positions=True)
		self.close_calculation()
		
	def molecular_dynamics(self, atoms, timestep, fix, step_iterator, update_cell, total_steps, constraints):
		fix = fix
		timestep = str(self.from_ase_units(timestep, 'time'))
		self.setup_calculation(atoms)
		self.evaluate_forces()
		cur_step = 0
		
		f = self.lammps_process
		f.write('fix mdfix %s\n' % fix)
		f.write('timestep %s\n' % timestep)
		for c in constraints:
			for cmd in c.get_commands(self.atoms):
				f.write(cmd + '\n')
		
		for nsteps in step_iterator:
			self.set_dumpfreq(cur_step+nsteps)
			if total_steps:
				f.write('run %s start 0 stop %s\n' % (nsteps, total_steps))
			else:
				f.write('run %s\n' % nsteps)
			f.write('print "%s"\n' % CALCULATION_END_MARK)
			f.flush()
			self.lammps_process.read_lammps_output()
			self.read_lammps_trj(update_positions=True, update_cell=update_cell)
			
			self.atoms_after_last_calc = atoms
			cur_step += nsteps
			
		self.close_calculation()
	
	def evaluate_forces(self):
		self.set_dumpfreq(1)
		f = self.lammps_process
		f.write('run 0\n')
		f.write('print "%s"\n' % CALCULATION_END_MARK)
		f.flush()
		self.lammps_process.read_lammps_output()
		self.read_lammps_trj(update_positions=False)
	
	def setup_calculation(self, atoms):
		filelabel = self.prepare_lammps_io()
		
		if not self.lammps_process.running():
			self.lammps_process.start(self.tmp_dir, self.lammps_command, filelabel)
				
		if np.all(atoms.pbc == False):
			# Make sure the atoms are inside the cell
			inv_cell = np.linalg.inv(atoms.cell)
			frac_positions = np.dot(inv_cell, atoms.positions.T)
			if np.any(frac_positions < 0) or np.any(frac_positions > 1):
				atoms.center(vacuum=0)

		self.atoms = atoms
		self.prism = Prism(self.atoms.cell)
		self.prepare_data()
		self.prepare_calculation(self.atoms, self.data)
		self.write_lammps_input()
		self.calls += 1
	
	def close_calculation(self):
		if not self.keep_alive:
			self.lammps_process.terminate()

		exitcode = self.lammps_process.poll()
		if exitcode and exitcode != 0:
			raise RuntimeError('LAMMPS exited in %s with exit code: %d.' %\
								(self.tmp_dir, exitcode))
		
		self.lammps_trj_file.close()
		self.lammps_inputdata_file.close()
		if self.debug == True: self.lammps_process.close_logs()
		
			
	def prepare_lammps_io(self):
		label = '%s%06d' % (self.label, self.calls)
		
		# In python 2.4 the delete kwarg doesn't exist
		#args = dict(dir=self.tmp_dir, delete=(not self.debug))
		#self.lammps_trj_file = NamedTemporaryFile(mode='r', prefix='trj_'+label, **args)
		#self.lammps_inputdata_file = NamedTemporaryFile(prefix='data_'+label, **args)
		
		if self.debug:
			self.lammps_trj_file = open(os.path.join(self.tmp_dir, 'trj_'+label), 'w')
			self.lammps_trj_file.close()
			self.lammps_trj_file = open(self.lammps_trj_file.name)
			self.lammps_inputdata_file = open(os.path.join(self.tmp_dir, 'data_'+label), 'w')
		else:
			self.lammps_trj_file = NamedTemporaryFile(mode='r', prefix='trj_'+label, dir=self.tmp_dir)
			self.lammps_inputdata_file = NamedTemporaryFile(prefix='data_'+label, dir=self.tmp_dir)
			
		return label
	
	
	def prepare_data(self):
		""" Prepare self.data for write_lammps_data() using self.ff_data """
		
		atoms = self.atoms
		tables = []
		ff_data = self.ff_data
		
		if not 'bonds' in atoms.info:
			atoms.info['bonds'] = Bonds(atoms, autodetect=True)
		
		atom_types = self.atom_types(atoms)
		atom_typeorder = list(set(atom_types))
		
		b = atoms.info['bonds']
		bonds = []
		if self.parameters.bond_style:
			for indices in b:
				type = SequenceType([atom_types[i] for i in indices])
				actual_type = ff_data.get_actual_type('bond', type)
				bonds.append(dict(indices=indices, type=actual_type))
			
		angles = []
		if self.parameters.angle_style:
			for indices in b.find_angles():
				type = SequenceType([atom_types[i] for i in indices])
				actual_type = ff_data.get_actual_type('angle', type)
				angles.append(dict(indices=indices, type=actual_type))
			
		dihedrals = []
		if self.parameters.dihedral_style:
			for indices in b.find_dihedrals():
				type = SequenceType([atom_types[i] for i in indices])
				actual_type = ff_data.get_actual_type('dihedral', type)
				dihedrals.append(dict(indices=indices, type=actual_type))
		
		impropers = []
		if self.parameters.improper_style:
			for indices in b.find_impropers():
				if ff_data.class2:
					j, i, k, l = indices
				else:
					i, j, k, l = indices
				a, c, b, d = (atom_types[ind] for ind in indices)
				type = ImproperType(central_type=a, other_types=(c,b,d), class2=ff_data.class2)
				actual_type = ff_data.get_actual_type('improper', type)
				impropers.append(dict(indices=(i,j,k,l), type=actual_type))
		
		# Coeffs	
		def add_coeff_tables(param_group, objects, typeorder=None, warn_missing=True):
			if not objects: return
			if typeorder:
				used_types = typeorder
			else:
				used_types = set(object['type'] for object in objects)
			
			typeorder = []
			available_tables = ff_data.available_tables(param_group)
			new_tables = {}
			for type in used_types:
				try:
					actual_type = ff_data.get_actual_type(param_group, type)
					params = ff_data.get(param_group, actual_type)
				except KeyError:
					if warn_missing: print 'No parameters for %s!' % actual_type
					continue
				typeorder.append(actual_type)
				
				for title, ncols in available_tables:
					try:
						values = params[title]
					except KeyError:
						if warn_missing: print 'No %s for %s!' % (title, actual_type)
						values = [0]*ncols
					table = new_tables.setdefault(title, [])
					table.append(values)
			tables.extend(new_tables.items())
			return typeorder
		
		# Add masses to ff_data
		masses = dict(zip(atom_types, self.atoms.get_masses()))
		for type in atom_typeorder:
			actual_type = ff_data.get_actual_type('atom', type)
			ff_data.add('atom', actual_type, 'Masses', [masses[type]])
			
		add_coeff_tables('atom', atom_types, atom_typeorder)
		
		bond_typeorder     = add_coeff_tables('bond', bonds, warn_missing=self.debug)
		angle_typeorder    = add_coeff_tables('angle', angles, warn_missing=self.debug)
		dihedral_typeorder = add_coeff_tables('dihedral', dihedrals, warn_missing=self.debug)
		improper_typeorder = add_coeff_tables('improper', impropers, warn_missing=self.debug)
		
		# Atoms
		self.set_charges(atoms, atom_types)
		atom_typeids = [atom_typeorder.index(at)+1 for at in atom_types]
		charges = self.atoms.get_charges()
		positions = self.prism.vector_to_lammps(self.atoms.positions)
		positions = self.from_ase_units(positions, 'distance')
		columns = [atom_typeids, charges, positions[:,0], positions[:,1], positions[:,2]]
		
		if self.parameters.atom_style == 'full':
			columns.insert(0, ['1']*len(self.atoms))
		elif not self.parameters.atom_style == 'charge':
			raise RuntimeError('Unsupported atom_style: %s' % self.parameters.atom_style)
			
		tables.append(('Atoms', zip(*columns)))
		
		# Bonds, Angles, etc.
		def add_object_table(title, objects, typeorder):
			if not objects or not typeorder: return
			table = []
			used_objects = []
			for obj in objects:
				if obj['type'] not in typeorder:
					continue
				typeid = typeorder.index(obj['type'])+1
				atoms = [idx+1 for idx in obj['indices']]
				table.append([typeid] + atoms)
				used_objects.append(obj)
			tables.append((title, table))
			return used_objects
		
		bonds     = add_object_table('Bonds', bonds, bond_typeorder)
		angles    = add_object_table('Angles', angles, angle_typeorder)
		dihedrals = add_object_table('Dihedrals', dihedrals, dihedral_typeorder)
		impropers = add_object_table('Impropers', impropers, improper_typeorder)
		
		if self.atoms.has('momenta'):
			vel = self.prism.vector_to_lammps(self.atoms.get_velocities())
			lammps_velocities = self.from_ase_units(vel, 'velocity')
			tables.append(('Velocities', lammps_velocities))
		
		data = self.data
		data.tables = tables
		data.atom_types = atom_types
		data.bonds      = bonds
		data.angles     = angles
		data.dihedrals  = dihedrals
		data.impropers  = impropers
		data.atom_typeorder     = atom_typeorder
		data.bond_typeorder     = bond_typeorder
		data.angle_typeorder    = angle_typeorder
		data.dihedral_typeorder = dihedral_typeorder
		data.improper_typeorder = improper_typeorder
		
	
	def write_lammps_input(self):
		"""Write LAMMPS parameters and input data """
		f = self.lammps_process
		parameters = self.parameters
		
		self.write_lammps_data()
		
		f.write('# (written by ASE)\n')
		f.write('clear\n')
		f.write('atom_style %s \n' % parameters.atom_style)
		
		pbc = self.atoms.get_pbc()
		f.write('units %s \n' % parameters.units)
		f.write('boundary %s %s %s \n' % tuple('sp'[int(x)] for x in pbc))
		if parameters.neighbor:
			f.write('neighbor %s \n' % (parameters.neighbor))
		if parameters.newton:
			f.write('newton %s \n' % (parameters.newton))

		# Write interaction stuff
		f.write('\n### interactions \n')
		if parameters.pair_style:
			f.write('pair_style %s \n' % parameters.pair_style)
			
		if parameters.bond_style and self.data.bonds:
			f.write('bond_style %s \n' % parameters.bond_style)
			
		if parameters.angle_style and self.data.angles:
			f.write('angle_style %s \n' % parameters.angle_style)
			
		if parameters.dihedral_style and self.data.dihedrals:
			f.write('dihedral_style %s \n' % parameters.dihedral_style)
			
		if parameters.improper_style and self.data.impropers:
			f.write('improper_style %s \n' % parameters.improper_style)
		
		if parameters.kspace_style:
			f.write('kspace_style %s \n' % parameters.kspace_style)
		
		if parameters.special_bonds:
			f.write('special_bonds %s \n' % parameters.special_bonds)
			
		if parameters.pair_modify:
			f.write('pair_modify %s \n' % parameters.pair_modify)
		
		f.write('\n### read data \n')
		f.write('read_data %s\n' % self.lammps_inputdata_file.name)
		
		# Extra pair coeffs
		for line in parameters.pair_coeffs:
			f.write('pair_coeff %s \n' % line)
		
		# Create groups
		for group_id, indices in parameters.groups:
			indices_str = ' '.join([str(i+1) for i in indices])
			f.write('group %s id %s\n' % (group_id, indices_str))
		
		for cmd in parameters.extra_cmds:
			f.write(cmd + '\n')
		
		f.write('\nthermo_style custom %s\n' % (' '.join(self._custom_thermo_args)))
		f.write('thermo 0\n')
		
		f.write('\ndump dump_all all custom ' +
				'1 %s %s\n' % (self.lammps_trj_file.name, ' '.join(self._dump_fields)) )

				
	def set_dumpfreq(self, freq):
		self.lammps_process.write('dump_modify dump_all every %s\n' % freq)
		
		
	def write_lammps_data(self):
		"""Write system configuration and force field parameters to file to be read
		with read_data by LAMMPS."""
		f = self.lammps_inputdata_file
		data = self.data
		
		f.write(f.name + ' (written by ASE) \n\n')

		f.write('%d \t atoms \n' % len(data.atom_types))
		if data.bonds:     f.write('%d \t bonds \n' % len(data.bonds))
		if data.angles:    f.write('%d \t angles \n' % len(data.angles))
		if data.dihedrals: f.write('%d \t dihedrals \n' % len(data.dihedrals))
		if data.impropers: f.write('%d \t impropers \n' % len(data.impropers))
		
		f.write('%d  atom types\n' % len(data.atom_typeorder))
		if data.bonds:     f.write('%d  bond types\n' % len(data.bond_typeorder))
		if data.angles:    f.write('%d  angle types\n' % len(data.angle_typeorder))
		if data.dihedrals: f.write('%d  dihedral types\n' % len(data.dihedral_typeorder))
		if data.impropers: f.write('%d  improper types\n' % len(data.improper_typeorder))

		xhi, yhi, zhi, xy, xz, yz = self.prism.get_lammps_prism()
		f.write('0.0 %f  xlo xhi\n' % xhi)
		f.write('0.0 %f  ylo yhi\n' % yhi)
		f.write('0.0 %f  zlo zhi\n' % zhi)
		
		if self.prism.is_skewed():
			f.write('%f %f %f  xy xz yz\n' % (xy, xz, yz))
		f.write('\n\n')
		
		for title, table in data.tables:
			if len(table) == 0:  continue
			f.write('%s \n\n' % title)
			for index, row in enumerate(table):
				f.write(('%d'+' %s'*len(row) +'\n') % ((index+1,) + tuple(row)))
			f.write('\n\n')
		
		f.flush()


	def read_lammps_trj(self, update_positions=False, update_cell=False):
		dump = read_lammps_dump(self.lammps_trj_file)
		
		rotate = self.prism.vector_to_ase
		self.forces = rotate(self.to_ase_units(dump.info['forces'], 'force'))
		
		if update_positions:
			dump.positions -= dump.info['celldisp']
			
			self.atoms.positions = rotate(self.to_ase_units(dump.positions, 'distance'))
			self.atoms.set_velocities(rotate(self.to_ase_units(dump.get_velocities(), 'velocity')))			
			if np.isnan(self.atoms.positions).any():
					raise RuntimeError('NaN detected in atomic coordinates!')
				
			if update_cell:
				self.atoms.set_cell(self.prism.vector_to_ase(dump.cell))
				
		self.atoms_after_last_calc = self.atoms.copy()
			
		if self.update_charges:
			self.atoms.set_charges(dump.get_charges())
	
	
	def to_ase_units(self, value, quantity):
		return unitconversion.convert(value, quantity, self.parameters.units, 'ASE')

	def from_ase_units(self, value, quantity):
		return unitconversion.convert(value, quantity, 'ASE', self.parameters.units)


class LammpsProcess:
	""" A class to handle the lammps process and read thermo output. There are
		sometimes errors related to the communication with the process and it
		is advisable to restart lammps after every calculation.
	"""
	def __init__(self, log=False):
		self.inlog = None
		self.outlog = None
		self.proc = None
		self.log = log
		
		
		self.thermo_output = []
		
	def __del__(self):
		if self.running(): self.proc.terminate()
		
	def invoke_lammps(self, tmp_dir, lammps_command, filelabel):
		if not lammps_command:
			lammps_command = os.environ.get('LAMMPS_COMMAND')
		if not lammps_command or len(lammps_command.strip()) == 0:
			raise RuntimeError('Please set LAMMPS_COMMAND environment variable')
		
		lammps_cmd_line = shlex.split(lammps_command)
		
		# Make sure we execute using the absolute path              
		lammps_cmd_line[0] = os.path.abspath(lammps_cmd_line[0])
		
		lammps_cmd_line += ['-log', '/dev/null']
		
		if self.log == True:
			# Save LAMMPS input and output for reference
			
			# in python 2.4 the delete=False option doesn't exist
			#self.inlog  = NamedTemporaryFile(prefix='in_'+filelabel, dir=tmp_dir, delete=False)
			#self.outlog = NamedTemporaryFile(prefix='log_'+filelabel, dir=tmp_dir, delete=False)
			self.inlog  = open(os.path.join(tmp_dir, 'in_'+filelabel), 'w')
			self.outlog = open(os.path.join(tmp_dir, 'out_'+filelabel), 'w')
			
		return Popen(lammps_cmd_line,
							cwd=tmp_dir, stdin=PIPE, stdout=PIPE, stderr=sys.stderr)
	
	def start(self, tmp_dir, lammps_command=None, filelabel=''):
		if self.running(): self.terminate()
		self.proc = self.invoke_lammps(tmp_dir, lammps_command, filelabel)
	
	def running(self):
		return self.proc and self.proc.poll() == None
	
	def poll(self):
		return self.proc.poll()
	
	def terminate(self):
		if not self.running(): pass
		self.proc.stdin.close()
		return self.proc.wait()
	
	def write(self, data):
		self.proc.stdin.write(data)
		if self.inlog: self.inlog.write(data)
	
	def readline(self):
		line = self.proc.stdout.readline()
		if self.outlog: self.outlog.write(line)
		return line
		
	def flush(self):
		self.proc.stdin.flush()
		self.write('log /dev/null\n')
		
	def close_logs(self):
		if self.inlog: self.inlog.close()
		if self.outlog: self.outlog.close()
		self.inlog = None
		self.outlog = None
	
	
	def get_thermo(self, key):
		""" Return the value of thermo variable key """
		return self.thermo_output[-1][key]
		
	def read_lammps_output(self):
		""" Read thermo output from LAMMPS stdout """
		f = self
		
		def translate_keys(keys):
			result = []
			for key in keys:
				if key == 'PotEng': k = 'pe'
				elif key == 'KinEng': k = 'ke'
				elif key == 'TotEng': k = 'etotal'
				else: k = key.lower()
				result.append(k)
			return result
		
		thermo_output = []
		line = f.readline()
		while line and line.strip() != CALCULATION_END_MARK:
			if 'ERROR' in line:
				raise RuntimeError('LAMMPS execution failed. LAMMPS %s' % line)
			
			words = line.split()
			if words and np.all([w[0].isupper() for w in words]):
				# Seems to be the start of thermo output
				keys = translate_keys(words)
				while True:
					line = f.readline()
					fields = line.split()
					if len(fields) != len(keys): break
					try:
						fields = map(float, fields) # convert to float
						thermo_output.append(dict(zip(keys, fields)))
					except ValueError:
						break   # Wasn't a thermo line after all
			else:
				line = f.readline()

		self.thermo_output = thermo_output
		if len(thermo_output) == 0:
			raise RuntimeError('No thermo output from LAMMPS!')


class Prism:
	"""The representation of the unit cell in LAMMPS"""

	def __init__(self, cell, pbc=(True,True,True), digits=10):
		# Use LQ decomposition to get the lammps cell
		# ase_cell * R = lammps_cell
		Qtrans, Ltrans = np.linalg.qr(cell.T)
		self.R = Qtrans
		self.lammps_cell = Ltrans.T
	
		if self.is_skewed() and not np.all(pbc):
			raise RuntimeError('Skewed lammps cells MUST have '
							'PBC == True in all directions!')

	def get_lammps_prism(self):
		return self.lammps_cell[(0,1,2,1,2,2), (0,1,2,0,0,1)]
		
	def update_cell(self, xyz, offdiag):
		self.lammps_cell = self.to_cell_matrix(xyz, offdiag)
		return np.dot(self.lammps_cell, self.R.T)
	
	def to_cell_matrix(self, xyz, offdiag):
		x, y, z = xyz
		xy, xz, yz = offdiag
		return np.array([[x,0,0], [xy,y,0], [xz, yz, z]])

	def vector_to_lammps(self, vec):
		return np.dot(vec, self.R)
		
	def vector_to_ase(self, vec):
		return np.dot(vec, self.R.T)
	
	def is_skewed(self):
		tolerance = 1e-6
		cell_sq = self.lammps_cell**2
		return np.sum(np.tril(cell_sq, -1)) / np.sum(np.diag(cell_sq)) > tolerance


