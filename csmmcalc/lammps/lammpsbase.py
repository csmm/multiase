#!/usr/bin/env python

# An ASE calculator for the LAMMPS classical MD code available from
#       http://lammps.sandia.gov/
# The environment variable LAMMPS_COMMAND must be defined to point to the LAMMPS binary.
#
# Copyright (C) 2012 Tuukka Verho, tuukka.verho@aalto.fi
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this file; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# or see <http://www.gnu.org/licenses/>.


import os, sys
import shutil
import shlex
import time
from subprocess import Popen, PIPE
from threading import Thread
from tempfile import mkdtemp, NamedTemporaryFile
import numpy as np
import decimal as dec
from ase import Atoms
import lammpsUnitConversion
import StringIO
import threading, Queue

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
		self.masses = None
		self.coeff_tables = []
		self.atom_types = None
		self.bonds = None
		self.angles = None
		self.dihedrals = None
		self.impropers = None
		
	def add_coeffs(self, title, table):
		self.coeff_tables.insert(dict(title=title, table=table))
		
		
		
class LAMMPSBase:

	def __init__(self, label='lammps', tmp_dir=None, parameters={}, 
				update_charges = False, lammps_command=None,
				keep_alive=False, debug=False):
		"""The LAMMPS calculators object """

		self.label = label
		self.parameters = LAMMPSParameters(**parameters)
		self.data = LAMMPSData()
		self.calls = 0
		self.forces = None
		self.atoms = None
		self.atoms_after_last_calc = None
		self.error = None
		self.update_charges = update_charges
		self.keep_alive = keep_alive
		self.lammps_command = lammps_command
		self.debug = debug           
		self.lammps_process = LammpsProcess()
		
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
	
	def prepare_calculation(self):
		""" Implement this method in subclasses """
		raise NotImplementedException()
		
		
	def use_ff_file(self, filepath, target_filename=None):
		if not target_filename: 
			target_filename=os.path.split(filepath)[1]
		shutil.copy(filepath, os.path.join(self.tmp_dir, target_filename))       

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
		self.prepare_calculation()
		self.prism = prism(self.atoms.get_cell())
		self.prepare_lammps_io()
		self.write_lammps_input()
	
	
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
		f.write('log /dev/stdout\n') # Force LAMMPS to flush log
		f.flush()
		
		
	def run_minimization(self, param_string, min_style=None):
		self.set_dumpfreq(999999)
		f = self.lammps_process
		if min_style:
			f.write('min_style %s\n' % min_style)
		f.write('minimize %s\n'  % param_string)
		f.write('print "%s"\n' % CALCULATION_END_MARK)
		f.write('log /dev/stdout\n') # Force LAMMPS to flush log
		
	def run_md(self, fix, timestep, nsteps):
		f = self.lammps_process
		f.write('fix fix_all %s\n' % fix)
		f.write('timestep %s\n' % timestep)
		f.write('run %s\n' % nsteps)
		f.write('print "%s"\n' % CALCULATION_END_MARK)
		f.write('log /dev/stdout\n') # Force LAMMPS to flush log
	
			
	def write_lammps_input(self):
		"""Write LAMMPS parameters and input data """
		f = self.lammps_process
		parameters = self.parameters
		
		f.write('# (written by ASE)\n')
		f.write('clear\n')
		f.write('atom_style %s \n' % parameters.atom_style)
		
		pbc = self.atoms.get_pbc()
		f.write('units %s \n' % parameters.units)
		f.write('boundary %c %c %c \n' % tuple('sp'[x] for x in pbc))
		if parameters.neighbor:
			f.write('neighbor %s \n' % (parameters.neighbor))
		if parameters.newton:
			f.write('newton %s \n' % (parameters.newton))

		# Write interaction stuff
		f.write('\n### interactions \n')
		if parameters.pair_style:
			f.write('pair_style %s \n' % parameters.pair_style)
			
		if parameters.bond_style and len(self.data.bonds) != 0:
			f.write('bond_style %s \n' % parameters.bond_style)
			
		if parameters.angle_style and len(self.data.angles) != 0:
			f.write('angle_style %s \n' % parameters.angle_style)
			
		if parameters.dihedral_style and len(self.data.dihedrals) != 0:
			f.write('dihedral_style %s \n' % parameters.dihedral_style)
			
		if parameters.improper_style and len(self.data.impropers) != 0:
			f.write('improper_style %s \n' % parameters.improper_style)
		
		if parameters.special_bonds:
			f.write('special_bonds %s \n' % parameters.special_bonds)
			
		if parameters.pair_modify:
			f.write('pair_modify %s \n' % parameters.pair_modify)
		
		self.write_lammps_data()
		f.write('\n### read data \n')
		f.write('read_data %s\n' % self.lammps_inputdata_file.name)
		
		# Extra pair coeffs
		for line in parameters.pair_coeffs:
			f.write('pair_coeff %s \n' % line)
		
		for cmd in parameters.extra_cmds:
			f.write(cmd + '\n')
		
		f.write(('thermo_style custom %s\n' +
				'thermo_modify flush yes\n' +
				'thermo 0\n') % (' '.join(self._custom_thermo_args)))
				
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
		
		def create_set_of_types(lst):
			if not lst: return None
			return set(element[0] for element in lst)
		
		atom_types     = set(self.data.atom_types)
		bond_types     = create_set_of_types(data.bonds)
		angle_types    = create_set_of_types(data.angles)
		dihedral_types = create_set_of_types(data.dihedrals)
		improper_types = create_set_of_types(data.impropers)

		f.write('%d  atom types\n' % len(atom_types))
		if data.bonds:     f.write('%d  bond types\n' % len(bond_types))
		if data.angles:    f.write('%d  angle types\n' % len(angle_types))
		if data.dihedrals: f.write('%d  dihedral types\n' % len(dihedral_types))
		if data.impropers: f.write('%d  improper types\n' % len(improper_types))

		p = prism(self.atoms.get_cell())
		xhi, yhi, zhi, xy, xz, yz = p.get_lammps_prism()

		f.write('0.0 %f  xlo xhi\n' % xhi)
		f.write('0.0 %f  ylo yhi\n' % yhi)
		f.write('0.0 %f  zlo zhi\n' % zhi)
		
		if p.is_skewed():
			f.write('%f %f %f  xy xz yz\n' % (xy, xz, yz))
		f.write('\n\n')
	
		def print_table(title, table):
			if len(table) == 0: return
			f.write('%s \n\n' % title)
			for index, row in enumerate(table, 1):
				f.write(('%d'+' %s'*len(row) +'\n') % ((index,) + tuple(row)))
			f.write('\n\n')
		
		if data.masses:
			print_table('Masses', [[mass] for mass in data.masses])
		
		for coeff_table in self.data.coeff_tables:
			print_table(coeff_table['title'], coeff_table['data'])
		
		# Atoms
		charges = self.atoms.get_charges()
		lammps_positions = p.pos_to_lammps(self.atoms.positions)
		positions = self.from_ase_units(lammps_positions, 'distance')
		if self.parameters.atom_style == 'full':
			table = [['1', tp, q] + list(pos)
				for tp, q, pos in zip(self.data.atom_types, charges, positions)]
		elif  self.parameters.atom_style == 'charge':
			table = [[tp, q] + list(pos)
				for tp, q, pos in zip(self.data.atom_types, charges, positions)]
		else:
			raise RuntimeError('Unsupported atom_style: %s' % self.parameters.atom_style)
		print_table('Atoms', table)
		
		if data.bonds:
			print_table('Bonds', data.bonds)
		if data.angles:
			print_table('Angles', self.data.angles)
		if data.dihedrals:
			print_table('Dihedrals', self.data.dihedrals)
		if data.impropers:
			print_table('Impropers', self.data.impropers)
		
		if self.atoms.has('momenta'):
			vel = self.atoms.get_velocities()
			lammps_velocities = self.from_ase_units(vel, 'velocity')
			print_table('Velocities', lammps_velocities)
		
		f.flush()
		
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
				if not line: print 'Line is', repr(line)
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
				lo = [] ; hi = [] ; tilt = []
				id = [] ; type = []
				positions  = np.empty((len(self.atoms), 3)) * np.nan
				forces     = np.empty((len(self.atoms), 3)) * np.nan
				velocities = np.empty((len(self.atoms), 3)) * np.nan
				charges    = np.empty((len(self.atoms), )) * np.nan

			if 'ITEM: NUMBER OF ATOMS' in line:
				line = f.readline()
				n_atoms = int(line.split()[0])
			
			if 'ITEM: BOX BOUNDS' in line:
				# save labels behind "ITEM: BOX BOUNDS" in triclinic case (>=lammps-7Jul09)
				tilt_items = line.split()[3:]
				for i in range(3):
					line = f.readline()
					fields = line.split()
					lo.append(float(fields[0]))
					hi.append(float(fields[1]))
					if (len(fields) >= 3):
						tilt.append(float(fields[2]))
			
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
										
					rotate = self.prism.pos_to_ase
					positions[id-1]  = rotate(self.to_ase_units(np.array(pos), 'distance'))
					velocities[id-1] = rotate(self.to_ase_units(np.array(vel), 'velocity'))
					forces[id-1]     = rotate(self.to_ase_units(np.array(force), 'force'))
					charges[id-1]    = self.to_ase_units(np.array(q), 'charge')
				
				self.forces = forces
				if timestep > 0:
					self.atoms.set_positions(positions)
					self.atoms.set_velocities(velocities)
					
					if len(tilt) == 0:
						tilt = [0, 0, 0]
					if (np.array(hi + tilt) != self.prism.get_lammps_prism()).any():
						# Update unit cell
						self.atoms.cell = self.prism.update_cell(hi, tilt)
					
				if self.update_charges:
					self.atoms.set_charges(charges)
				
				if np.isnan(positions).any():
					raise RuntimeError('NaN detected in atomic coordinates!')
	
	
	def to_ase_units(self, value, quantity):
		return lammpsUnitConversion.convert(value, quantity, self.parameters.units, 'ASE')

	def from_ase_units(self, value, quantity):
		return lammpsUnitConversion.convert(value, quantity, 'ASE', self.parameters.units)


class LammpsProcess:
	""" A class to handle the lammps process. There are sometimes errors related
	    to the communication with the process and it is advisable to restart lammps
	    after every calculation.
	"""
	def __init__(self, inlog=None, outlog=None):
		self.inlog = inlog
		self.outlog = outlog
		self.proc = None
		
	def __del__(self):
		if self.running(): self.proc.terminate()
		
	def invoke_lammps(self, tmp_dir, lammps_command):
		if not lammps_command:
			lammps_command = os.environ.get('LAMMPS_COMMAND')
		if not lammps_command or len(lammps_command.strip()) == 0:
			raise RuntimeError('Please set LAMMPS_COMMAND environment variable')
		
		lammps_cmd_line = shlex.split(os.environ['LAMMPS_COMMAND'])
		
		# Make sure we execute using the absolute path              
		lammps_cmd_line[0] = os.path.abspath(lammps_cmd_line[0])
		
		stdout = NamedTemporaryFile(mode='w')
		self.tmp = stdout
		
		self.stdout = open(stdout.name, 'r')
		
		# Because LAMMPS output to stdout does not work with subprocess
		# (for unknown reason), we direct logging output to stdout as a workaround.
		return Popen(lammps_cmd_line +['-log', '/dev/stdout'],
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
		
	def update_cell(self, hi, tilt=(0,0,0)):
		xhi, yhi, zhi = hi
		xy, xz, yz = tilt
		self.lammps_cell = np.array([[xhi,0,0], [xy,yhi,0], [xz, yz, zhi]])
		return self.pos_to_ase(self.lammps_cell)

	def pos_to_lammps(self, position):
		return np.dot(position, self.R)
	
	def pos_to_ase(self, position):
		return np.dot(position, self.R.T)
	
	def is_skewed(self):
		tolerance = 1e-6
		cell_sq = self.lammps_cell**2
		return np.sum(np.tril(cell_sq, -1)) / np.sum(np.diag(cell_sq)) > tolerance
		