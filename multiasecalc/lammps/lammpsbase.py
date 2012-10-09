import os, sys
import shutil, shlex
from subprocess import Popen, PIPE
from tempfile import mkdtemp, NamedTemporaryFile
import numpy as np
import unitconversion
from itertools import combinations, permutations

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
		self.special_bonds  = pars.get('special_bonds')
		self.pair_modify    = pars.get('pair_modify')
		
		# Pair coeffs to be specified in the input file
		self.pair_coeffs    = pars.get('pair_coeffs', [])
		
		self.extra_cmds     = pars.get('extra_cmds', [])

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

class FFData:
	def __init__(self):
		self.atom = {}
		self.bond = {}
		self.angle = {}
		self.dihedral = {}
		self.improper = {}
		
	def add(self, group, type, title, values):
		groupdict = getattr(self, group)
		d = groupdict.setdefault(type, {})
		d[title] = values


class LAMMPSBase:

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
		self.lammps_process = LammpsProcess()
		self.calls = 0
		
		# Thermo output of LAMMPS, a list of dictionaries. See read_lammps_output()
		self.thermo_content = []
		
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
		energy = self.thermo_content[-1]['pe']
		return self.to_ase_units(energy, 'energy')

	def get_forces(self, atoms):
		self.update(atoms)
		
		return self.forces

	def get_stress(self, atoms):
		self.update(atoms)
		tc = self.thermo_content[-1]
		stress =  np.array([tc[i] for i in ('pxx','pyy','pzz',
										'pyz','pxz','pxy')])
		return self.to_ase_units(stress, 'stress')  
	
	def calculation_required(self, atoms, quantities=None):
		return atoms != self.atoms_after_last_calc
		
	def update(self, atoms):
		if not self.calculation_required(atoms): return
		self.setup_calculation(atoms)
		self.run_single_step()
		self.read_lammps_output() 
		self.read_lammps_trj()
		self.close_calculation()
	
	def minimize(self, atoms, etol=0, ftol=0, maxeval=100000, maxiter=100000000, min_style=None):
		if etol == 0 and ftol == 0: 
			raise RuntimeError('Specify at least one tolerance value!')
		ftol = self.from_ase_units(ftol, 'force')
		minimize_params = '%g %g %i %i' % (etol, ftol, maxeval, maxiter)
		
		self.setup_calculation(atoms)
		
		# Observers may need data before the first iteration is run
		self.run_single_step()
		self.read_lammps_output()
		self.read_lammps_trj()
		
		self.run_minimization(minimize_params, min_style)
		self.read_lammps_output()
		self.read_lammps_trj()
		self.close_calculation()

		
	def molecular_dynamics(self, atoms, timestep, fix):
		nsteps = yield
		fix = fix
		timestep = str(self.from_ase_units(timestep, 'time'))
		self.setup_calculation(atoms)
		
		# Observers may need data before the first timestep is run
		self.run_single_step()
		self.read_lammps_output()
		self.read_lammps_trj()
		
		cur_step = 0
		while nsteps:
			self.set_dumpfreq(cur_step+nsteps)
			cur_step += nsteps
			self.run_md(fix, timestep, nsteps)
			self.read_lammps_output()
			self.read_lammps_trj()
			self.atoms_after_last_calc = atoms
			nsteps = yield
		self.close_calculation()
		
	def setup_calculation(self, atoms):
		if not self.lammps_process.running():
			self.lammps_process.start(self.tmp_dir, self.lammps_command)
				
		if all(atoms.pbc == False): 
			# Make sure the atoms are inside the cell
			inv_cell = np.linalg.inv(atoms.cell)
			frac_positions = np.dot(inv_cell, atoms.positions.T)
			if np.any(frac_positions < 0) or np.any(frac_positions > 1):
				atoms.center(vacuum=0)

		self.atoms = atoms
		if atoms != self.atoms_after_last_calc:
			self.prism = prism(self.atoms.get_cell())
			self.data.clear()
			self.prepare_data()
			self.prepare_calculation(self.atoms, self.data)
		self.prepare_lammps_io()
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
		
		self.atoms_after_last_calc = self.atoms.copy()
			
	def prepare_lammps_io(self):
		# setup input and output files for LAMMPS
		label = '%s%06d' % (self.label, self.calls)
	
		self.lammps_trj_file = NamedTemporaryFile(mode='r', prefix='trj_'+label,
									dir=self.tmp_dir, delete=(not self.debug))
		
		self.lammps_inputdata_file = NamedTemporaryFile(prefix='data_'+label,
									dir=self.tmp_dir, delete=(not self.debug))
	
		if self.debug == True:
			# Save LAMMPS input and output for reference
			inlog = NamedTemporaryFile(prefix='in_'+label, dir=self.tmp_dir, delete=False)
			outlog = NamedTemporaryFile(prefix='log_'+label, dir=self.tmp_dir, delete=False)
			self.lammps_process.inlog = inlog
			self.lammps_process.outlog = outlog
		elif self.debug == 'stdout':
			inlog = None
			outlog = sys.stdout
		
	def run_single_step(self):
		self.set_dumpfreq(1)
		f = self.lammps_process
		f.write('run 0\n')
		f.write('print "%s"\n' % CALCULATION_END_MARK)
		f.flush()
		
	def run_minimization(self, param_string, min_style=None):
		self.set_dumpfreq(999999)
		f = self.lammps_process
		if min_style:
			f.write('min_style %s\n' % min_style)
		f.write('minimize %s\n'  % param_string)
		f.write('print "%s"\n' % CALCULATION_END_MARK)
		f.flush()
		
	def run_md(self, fix, timestep, nsteps):
		f = self.lammps_process
		f.write('fix fix_all %s\n' % fix)
		f.write('timestep %s\n' % timestep)
		f.write('run %s\n' % nsteps)
		f.write('print "%s"\n' % CALCULATION_END_MARK)
		f.flush()
	
			
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
		f.write('boundary %s %s %s \n' % tuple('sp'[x] for x in pbc))
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
		
		if parameters.special_bonds:
			f.write('special_bonds %s \n' % parameters.special_bonds)
			
		if parameters.pair_modify:
			f.write('pair_modify %s \n' % parameters.pair_modify)
		
		f.write('\n### read data \n')
		f.write('read_data %s\n' % self.lammps_inputdata_file.name)
		
		# Extra pair coeffs
		for line in parameters.pair_coeffs:
			f.write('pair_coeff %s \n' % line)
		
		for cmd in parameters.extra_cmds:
			f.write(cmd + '\n')
		
		f.write('thermo_style custom %s\n' % (' '.join(self._custom_thermo_args)))
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
			#if len(table) == 0: return
			f.write('%s \n\n' % title)
			for index, row in enumerate(table, 1):
				f.write(('%d'+' %s'*len(row) +'\n') % ((index,) + tuple(row)))
			f.write('\n\n')
		
		f.flush()
	
	def prepare_data(self):
		atoms = self.atoms
		
		if not atoms.has('bonds'):
			detect_bonds(atoms)
			
		atom_types = self.atom_types(atoms)
		atom_typeorder = list(set(atom_types))
		bonds = Bond.all_bonds(atoms)
		angles = Angle.all_angles(atoms)
		dihedrals = Dihedral.all_dihedrals(atoms)
		impropers = Improper.all_impropers(atoms)
		
		for bond in bonds:
			bond.type = SequenceType(atom_types[i] for i in bond.atoms)
			
		for angle in angles:
			angle.type = SequenceType(atom_types[i] for i in angle.atoms)
			
		for dihedral in dihedrals:
			dihedral.type = SequenceType(atom_types[i] for i in dihedral.atoms)
			
		for improper in impropers:
			central_type = atom_types[improper.central_atom]
			other_types = [atom_types[i] for i in improper.other_atoms]
			improper.type = ImproperType([central_type] + other_types)
		
		tables = []
		ff_data = self.ff_data
		
		# Coeffs
		# Helper functions
		def get_tablenames(params):
			d = {}
			for entry in params.values():
				d.update((name, len(params)) for name, params in entry.items())
			return d.items()
			
			
		def coeff_table_generator(title, params, typeorder, empty_value, warn_missing):
			for type in typeorder:
				try:
					yield params[type][title]
				except TypeError:
					if warn_missing: print 'No %s parameters for %s!' % (title, type)
					yield empty_value
		
		def add_coeff_tables(params, objects, typeorder=None, warn_missing=False):
			if not params or not objects: return
			if not typeorder:
				used_types = set(object.type for object in objects)
				used_params = {}
				for type in used_types:
					try: used_params[type] = params[type]
					except KeyError:
						if warn_missing: print 'No parameters for %s!' % type
				params = used_params
				typeorder = params.keys()
				
			for title, ncols in get_tablenames(params):
				table = coeff_table_generator(title, params, typeorder, [0]*ncols, warn_missing)
				tables.append((title, table))
			return typeorder
		
		masses = dict(zip(atom_types, self.atoms.get_masses()))
		for type in atom_typeorder:
			ff_data.add('atom', type, 'Masses', [masses[type]])
			
		add_coeff_tables(ff_data.atom, atom_types, atom_typeorder)
		
		bond_typeorder     = add_coeff_tables(ff_data.bond, bonds, warn_missing=True)
		angle_typeorder    = add_coeff_tables(ff_data.angle, angles)
		dihedral_typeorder = add_coeff_tables(ff_data.dihedral, dihedrals)
		improper_typeorder = add_coeff_tables(ff_data.improper, impropers)
		
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
				if obj.type not in typeorder:
					objects.remove(obj)
					continue
				typeid = typeorder.index(obj.type)+1
				atoms = [idx+1 for idx in obj.atoms]
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
	
	
	def read_lammps_output(self):
		""" Read thermo output from LAMMPS stdout """
		f = self.lammps_process
		
		thermo_mark = ' '.join([x.capitalize() for x in
											self._custom_thermo_args[0:3]])
		
		thermo_content = []
		line = f.readline()
		while line and line.strip() != CALCULATION_END_MARK:
			if 'ERROR:' in line:
				raise RuntimeError('LAMMPS execution in %s failed. LAMMPS %s' 
						% (self.tmp_dir,line))
			if line.startswith(thermo_mark):
				# get thermo output
				while True:
					line = f.readline()
					fields = line.split()
					if len(fields) != len(self._custom_thermo_args): break
					try:
						fields = map(float, fields) # convert to float
					except ValueError:                       
						break   # Wasn't a thermo line after all
					thermo_content.append(dict(zip(self._custom_thermo_args, fields)))
					
			else:
				line = f.readline()
		self.thermo_content = thermo_content
		if len(thermo_content) == 0:
			raise RuntimeError('No thermo output from LAMMPS!')

		
	def read_lammps_trj(self):
		"""Read the LAMMPS dump file."""
		f = self.lammps_trj_file
		
		while True:
			line = f.readline()
			if not line: break

			if 'ITEM: TIMESTEP' in line:
				timestep = int(f.readline())
				n_atoms = 0
				lo_bound = [] ; hi_bound = [] ; tilt = []
				id = [] ; type = []
				positions  = np.ones((len(self.atoms), 3)) * np.nan
				forces     = np.ones((len(self.atoms), 3)) * np.nan
				velocities = np.ones((len(self.atoms), 3)) * np.nan
				charges    = np.ones((len(self.atoms), )) * np.nan

			if 'ITEM: NUMBER OF ATOMS' in line:
				line = f.readline()
				n_atoms = int(line.split()[0])
			
			if 'ITEM: BOX BOUNDS' in line:
				# save labels behind "ITEM: BOX BOUNDS" in triclinic case (>=lammps-7Jul09)
				tilt_items = line.split()[3:]
				for i in range(3):
					line = f.readline()
					fields = line.split()
					lo_bound.append(float(fields[0]))
					hi_bound.append(float(fields[1]))
					if (len(fields) >= 3):
						tilt.append(float(fields[2]))
						
				if len(tilt) == 0: 
					lo = np.array(lo_bound)
					hi = np.array(hi_bound)
					offdiag = [0,0,0]
				else:
					xy, xz, yz = tilt
					lo = np.array(lo_bound) - np.array([min(0,xy,xz,xy+xz), min(0,yz), 0])
					hi = np.array(hi_bound) - np.array([max(0,xy,xz,xy+xz), max(0,yz), 0])
					offdiag = tilt
			
			if 'ITEM: ATOMS' in line:
				# (reliably) identify values by labels behind "ITEM: ATOMS" - requires >=lammps-7Jul09
				column_names = line.split()[2:]
				for n in range(n_atoms):
					line = f.readline()
					fields = dict(zip(column_names, line.split()))
					id = int(fields['id'])
						
					pos   = [ float(fields[col]) for col in ['x', 'y', 'z'] ]
					vel   = [ float(fields[col]) for col in ['vx', 'vy', 'vz'] ]
					force = [ float(fields[col]) for col in ['fx', 'fy', 'fz'] ]
					q     = float(fields['q'])
					
					origin = lo
					rotate = self.prism.vector_to_ase
					positions[id-1]  = rotate(self.to_ase_units(np.array(pos)-origin, 'distance'))
					velocities[id-1] = rotate(self.to_ase_units(np.array(vel), 'velocity'))
					forces[id-1]     = rotate(self.to_ase_units(np.array(force), 'force'))
					charges[id-1]    = self.to_ase_units(q, 'charge')
				
				self.forces = forces
				if timestep > 0:
					self.atoms.set_positions(positions)
					self.atoms.set_velocities(velocities)
					
					self.atoms.cell = self.prism.update_cell(hi-lo, offdiag)
					
				if self.update_charges:
					self.atoms.set_charges(charges)
				
				if np.isnan(positions).any():
					raise RuntimeError('NaN detected in atomic coordinates!')
	
	
	def to_ase_units(self, value, quantity):
		return unitconversion.convert(value, quantity, self.parameters.units, 'ASE')

	def from_ase_units(self, value, quantity):
		return unitconversion.convert(value, quantity, 'ASE', self.parameters.units)


class LammpsProcess:
	""" A class to handle the lammps process. There are sometimes errors related
	    to the communication with the process and it is advisable to restart lammps
	    after every calculation.
	"""
	def __init__(self, inlog=None, outlog=None):
		self.inlog = inlog
		self.outlog = outlog
		self.proc = None
		self.output_hack = False  # see invoke_lammps()
		
	def __del__(self):
		if self.running(): self.proc.terminate()
		
	def invoke_lammps(self, tmp_dir, lammps_command):
		if not lammps_command:
			lammps_command = os.environ.get('LAMMPS_COMMAND')
		if not lammps_command or len(lammps_command.strip()) == 0:
			raise RuntimeError('Please set LAMMPS_COMMAND environment variable')
		
		lammps_cmd_line = shlex.split(lammps_command)
		
		# Make sure we execute using the absolute path              
		lammps_cmd_line[0] = os.path.abspath(lammps_cmd_line[0])
		
		if not 'mpirun' in lammps_cmd_line[0]:
			# If one doesn't execute LAMMPS with 'mpirun', the normal output to stdout
			# does not work. Here's a workaround.
			lammps_cmd_line += ['-log', '/dev/stdout']
			self.output_hack=True
		
		return Popen(lammps_cmd_line,
							cwd=tmp_dir, stdin=PIPE, stdout=PIPE, stderr=sys.stderr)
	
	def start(self, tmp_dir, lammps_command=None):
		if self.running(): self.terminate()
		self.proc = self.invoke_lammps(tmp_dir, lammps_command)
	
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
		if self.output_hack: self.write('log /dev/stdout\n')
		
	def close_logs(self):
		if self.inlog: self.inlog.close()
		if self.outlog: self.outlog.close()
		self.inlog = None
		self.outlog = None
	


class prism:
	"""The representation of the unit cell in LAMMPS"""

	def __init__(self, cell, pbc=(True,True,True), digits=10):
		# Use LQ decomposition to get the lammps cell
		# ase_cell * R = lammps_cell
		Qtrans, Ltrans = np.linalg.qr(cell.T)
		self.R = Qtrans
		self.lammps_cell = Ltrans.T
	
		if self.is_skewed() and not all(pbc):
			raise RuntimeError('Skewed lammps cells MUST have '
							'PBC == True in all directions!')

	def get_lammps_prism(self):
		return self.lammps_cell[(0,1,2,1,2,2), (0,1,2,0,0,1)]
	
	def cell_changed(self, xyz, offdiag):
		return (self.to_cell_matrix(xyz, offidiag) != self.lammps_cell).any()
	
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

		
def detect_bonds(atoms):
	from ase.data import covalent_radii
	tolerance = 1.1
	def bonded(i, j):
		bondLength = covalent_radii[atoms[i].number] + covalent_radii[atoms[j].number]
		return atoms.get_distance(i, j, mic=True) < bondLength*tolerance
			
	N = len(atoms)
	bonds = np.empty((N), dtype=object)
	for i in range(N):
		bonds[i] = [j for j in range(N) if bonded(i, j) and i != j]
	atoms.set_array('bonds', bonds)


class _Base:
	def __init__(self, atoms):
		self.atoms = list(atoms)
		self.type = None
	def __repr__(self):
		return str(self.type)
		
class Bond(_Base):
	@classmethod
	def all_bonds(cls, atoms):
		bondarr = atoms.get_array('bonds')
		return [Bond((a.index, b.index)) for a,b in combinations(atoms, 2) if b.index in bondarr[a.index]]

class Angle(_Base):
	@classmethod
	def all_angles(cls, atoms):
		angles = []
		for index, neigbors in enumerate(atoms.get_array('bonds')):
			for pair in combinations(neigbors, 2):
				angles.append(Angle((pair[0], index, pair[1])))
		return angles
	
class Dihedral(_Base):
	@classmethod
	def all_dihedrals(cls, atoms):
		dihedrals = []
		bonds = atoms.get_array('bonds')
		for j in range(len(atoms)):
			for i, k in permutations(bonds[j], 2):
				if k < j: continue
				for l in bonds[k]:
					if l == j: continue
					dihedrals.append(Dihedral((i,l,k,l)))
		return dihedrals
	
	
class Improper:
	def __init__(self, central_atom, other_atoms, class2=False):
		self.central_atom = central_atom
		self.other_atoms = list(other_atoms)
		self.type = None
		self.class2 = class2
	@property
	def atoms(self):
		if self.class2:
			return [other_atoms[0], self.central_atoms] + other_atoms[1:]
		else:
			return [self.central_atom] + other_atoms
			
	def __repr__(self):
		return '%s, %s' % (str(self.central_atom), str(self.other_atoms))
	@classmethod
	def all_impropers(cls, atoms):
		impropers = []
		for index, neighbors in enumerate(atoms.get_array('bonds')):
			for triplet in combinations(neighbors, 3):
				impropers.append(Improper(index, set(triplet)))
		return impropers
			 

class SequenceType:
	def __init__(self, atom_types):
		self.atom_types = sorted(atom_types)
	def __eq__(self, other):
		return self.atom_types == other.atom_types
	def __hash__(self):
		return sum(hash(tp) for tp in self.atom_types)
	def __repr__(self):
		return repr(self.atom_types)

class ImproperType:
	def __init__(self, atom_types, class2=False):
		if class2:
			self.central =  atom_types[0]
			self.others = sorted(atom_types[1:])
		else:
			self.central =  atom_types[1]
			self.others = sorted([atom_types[0]]+atom_types[2:])
	def get_types(self):
		return self.central, self.others
	def __eq__(self, other):
		return self.central == other.central and self.others == other.others
	def __hash__(self):
		return hash(self.central) + sum(hash(tp) for tp in self.others)
	def __repr__(self):
		return '%s, %s' %(self.central, self.others)
