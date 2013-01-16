from multiasecalc.lammps import unitconversion
from ase.optimize.optimize import Dynamics
from ase.io.trajectory import PickleTrajectory
from ase.md.logger import MDLogger
from ase import units
from random import random
import numpy as np

class LAMMPSOptimizer(Dynamics):
	""" Geometry optimizer for LAMMPS. works only with LAMMPS calculators """
	def __init__(self, atoms, restart=None, logfile=None, trajectory=None, algorithm='cg', relax_cell=False):
		Dynamics.__init__(self, atoms, logfile, trajectory)
		self.algorithm = algorithm
		self.relax_cell = relax_cell
		
	def run(self, fmax=0.001, steps=1e8):
		self.atoms.calc.minimize(self.atoms, ftol=fmax, maxeval=steps, min_style=self.algorithm, relax_cell=self.relax_cell)

		
class LAMMPSMolecularDynamics(Dynamics):
	""" Base class for molecular dynamics with LAMMPS. Requires a LAMMPS calculator. """
	def __init__(self, atoms, timestep, integrator='verlet', trajectory=None,
			traj_interval=1000, logfile=None, loginterval=100):
		Dynamics.__init__(self, atoms, None, None)

		self.dt = timestep
		
		if integrator == 'verlet':
			self.run_style = 'verlet'
		else:
			raise RuntimeError('Unknown integrator: %s' % thermostat)
		
		if trajectory:
			if isinstance(trajectory, str):
				trajectory = PickleTrajectory(trajectory, 'w', atoms)
			self.attach(trajectory, interval=traj_interval)
		
		if logfile:
			self.attach(MDLogger(dyn=self, atoms=atoms, logfile=logfile),
						interval=loginterval)
		
		self.fix = None
		self.cell_relaxed = False
		
	def run(self, steps=50, constraints=[]):
		self.nsteps = 0
		fix = 'all '+self.fix
		calc = self.atoms.calc
		it = self.run_iterator(steps)
		calc.molecular_dynamics(self.atoms, self.dt, fix, it, self.cell_relaxed, steps, constraints)
		
	def run_iterator(self, steps):
		cur_step = 0
		for target_step in range(steps+1):
			for function, interval, args, kwargs in self.observers:
				if target_step % interval == 0:
					if target_step > cur_step:
						yield target_step - cur_step
						cur_step = target_step
					function(*args, **kwargs)
		
		if cur_step < steps:
			yield steps - cur_step
    
	def get_time(self):
		return self.nsteps * self.dt

class LAMMPS_NVE(LAMMPSMolecularDynamics):
	""" Microcanonical ensemble """
	def __init__(self, atoms, timestep, **kwargs):
		LAMMPSMolecularDynamics.__init__(self, atoms, timestep, **kwargs)
		
		self.fix = 'nve'
		
class LAMMPS_NVT(LAMMPSMolecularDynamics):
	""" Constant temperature calculations with Nose-Hoover or Berendsen """
	def __init__(self, atoms, timestep, temperature, t_damp=100*units.fs,
				thermostat='Nose-Hoover', ramp_to_temp=None, **kwargs):
		LAMMPSMolecularDynamics.__init__(self, atoms, timestep, **kwargs)
		
		if thermostat == 'Nose-Hoover':
			cmd = 'nvt temp'
		elif thermostat == 'Berendsen':
			cmd = 'temp/berendsen'
		else:
			raise RuntimeError('Unknown thermostat: %s' % thermostat)
		
		t_damp = atoms.calc.from_ase_units(t_damp, 'time')
		
		if not ramp_to_temp: ramp_to_temp = temperature
		
		self.fix = '%s %f %f %f' %(cmd, temperature, ramp_to_temp, t_damp)
		

class LAMMPS_NPT(LAMMPSMolecularDynamics):
	""" Constant temperature and pressure calculations with Nose-Hoover """
	def __init__(self, atoms, timestep, temperature, externalstress, isotropic=True, t_damp=100*units.fs, p_damp=1000*units.fs, ramp_to_temp=None, **kwargs):
		LAMMPSMolecularDynamics.__init__(self, atoms, timestep, **kwargs)
			
		pressure = atoms.calc.from_ase_units(externalstress, 'pressure')
		t_damp = atoms.calc.from_ase_units(t_damp, 'time')
		p_damp = atoms.calc.from_ase_units(p_damp, 'time')
		
		if not ramp_to_temp: ramp_to_temp = temperature
		
		if hasattr(pressure, 'shape'):
			px, pxy, pxz = pressure[0,:]
			py, pyz = pressure[1,1:]
			pz = pressure[2,2]
			p_diags = [px, py, pz]
			args = ' '.join(['%s %f %f %f' % ('xyz'[i], p_diags[i], p_diags[i], p_damp) for i in range(3) if atoms.pbc[i]])
			if atoms.pbc[0] and atoms.pbc[1]:
				args += ' xy %f %f %f' % (pxy, pxy, p_damp)
			if atoms.pbc[1] and atoms.pbc[2]:
				args += ' yz %f %f %f' % (pyz, pyz, p_damp)
			if atoms.pbc[1] and atoms.pbc[2]:
				args += ' xz %f %f %f' % (pxz, pxz, p_damp)
		else:
			pvalues = '%f %f %f' % (pressure, pressure, p_damp)
		
			if atoms.pbc.all():
				if isotropic:
					coupling = 'iso'
				elif (np.dot(atoms.cell, atoms.cell) == atoms.cell**2).all():
					# orthogonal cell
					coupling = 'aniso'
				else:
					coupling = 'tri'
				args = '%s %s' % (coupling, pvalues)
			else:
				args = ' '.join(['%s %s' % ('xyz'[i], pvalues) for i in range(3) if atoms.pbc[i]])
				
		self.fix = 'npt temp %f %f %f %s' %(temperature, ramp_to_temp, t_damp, args)
		self.cell_relaxed = True

class SimpleConstraint:
	def __init__(self, indices):
		self.indices = indices
		
	def get_commands(self, atoms):
		fix = self.get_fix(atoms)
		id = '%s%s' % (self.__class__.__name__, abs(hash(tuple(self.indices))))
		groupname = 'group%s' % id
		fixname = 'fix%s' % id
		
		cmds = []
		indices_str = ' '.join([str(i+1) for i in self.indices])
		cmds.append('group %s id %s' % (groupname, indices_str))
		cmds.append('fix %s %s %s' % (fixname, groupname, fix))
		return cmds
		
	def get_fix(self, atoms):
		raise NotImplementedError()
		
class Spring(SimpleConstraint):
	def __init__(self, indices, point, spring_constant, R0=0.0):
		SimpleConstraint.__init__(self, indices)
		self.point = point
		self.K = spring_constant
		self.R0 = R0
		
	def get_fix(self, atoms):
		K = atoms.calc.from_ase_units(self.K, 'force')
		x, y, z = atoms.calc.prism.vector_to_lammps(self.point)
		return 'spring tether %f %f %f %f %f' % (K, x, y, z, self.R0)

class AddForce(SimpleConstraint):
	def __init__(self, indices, total_force):
		SimpleConstraint.__init__(self, indices)
		self.total_force = total_force
	
	def get_fix(self, atoms):
		force = self.total_force / len(self.indices)
		force = atoms.calc.prism.vector_to_lammps(force)
		fx, fy, fz = atoms.calc.from_ase_units(force, 'force')
		return 'addforce %f %f %f' % (fx, fy, fz)
		
		
class LJWall:
	def __init__(self, face, epsilon, sigma, wall_offset=None, final_wall_offset=None, mixing='arithmetic'):
		self.face = face
		self.epsilon = epsilon
		self.sigma = sigma
		self.offset = wall_offset
		self.final_offset = final_wall_offset
		self.mixing = mixing
		self.commands = []
		self.ngroups = 0
		self.nfixes = 0
		self.id = '%s%s' % (self.__class__.__name__, abs(hash(epsilon + sigma) + hash(face)))
		
		#if 'hi' in face:
		#	self.offset = -abs(self.offset)
		
	def get_commands(self, atoms):
		ffdata = atoms.calc.ff_data
		
		if self.final_offset != None:
			rampname = 'ramp%s' % self.id
			self.commands.append('variable %s equal ramp(%f,%f)' % (rampname, self.offset, self.final_offset))
			coord = 'v_%s' % rampname
		elif self.offset != None:
			coord = '%f' % self.offset
		else:
			coord = 'EDGE'
		
		for tp in atoms.calc.data.atom_typeorder:
			actual_type = ffdata.get_actual_type('atom', tp)
			eps, sig = ffdata.get_params('atom', actual_type)['Pair Coeffs']
			mixeps = np.sqrt(self.epsilon*eps)
			if self.mixing == 'arithmetic':
				mixsig = (self.sigma+sig)/2
			elif self.mixing == 'geometric':
				mixsig = np.sqrt(self.sigma*sig)
			else:
				raise RuntimeError('Invalid mixing type: %s' % self.mixing)
			typeid = atoms.calc.data.atom_typeorder.index(tp) + 1
			groupname = self.create_group_by_type(typeid)
			cutoff = 10.0
			fixstr = 'wall/lj126 %s %s %f %f %f units box pbc yes' % (self.face, coord, mixeps, mixsig, cutoff)
			self.create_fix(groupname, fixstr)
		
		return self.commands

	def create_group_by_type(self, typeid):
		groupname = 'group%s%s' % (self.id, typeid)
		self.commands.append('group %s type %i' % (groupname, typeid))
		self.ngroups += 1
		return groupname
		
	def create_fix(self, groupname, fixstr):
		fixname = 'fix%s%i' % (self.id, self.nfixes)
		self.commands.append('fix %s %s %s' % (fixname, groupname, fixstr))
		self.nfixes += 1
		return fixname
		
		