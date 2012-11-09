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
			traj_interval=1000, logfile=None, loginterval=100, constraints=[]):
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
		self.constraints = constraints
		
	def run(self, steps=50):
		self.nsteps = 0
		fix = 'all '+self.fix
		calc = self.atoms.calc
		it = self.run_iterator(steps)
		calc.molecular_dynamics(self.atoms, self.dt, fix, it, self.cell_relaxed, steps, self.constraints)
		
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
	def __init__(self, atoms, timestep, temperature, externalstress, t_damp=100*units.fs, p_damp=1000*units.fs, ramp_to_temp=None, **kwargs):
		LAMMPSMolecularDynamics.__init__(self, atoms, timestep, **kwargs)
			
		pressure = atoms.calc.from_ase_units(externalstress, 'pressure')
		t_damp = atoms.calc.from_ase_units(t_damp, 'time')
		p_damp = atoms.calc.from_ase_units(p_damp, 'time')
		
		if not ramp_to_temp: ramp_to_temp = temperature
		
		if (np.dot(atoms.cell, atoms.cell) == atoms.cell**2).all():
			# orthogonal
			coupling = 'aniso'
		else:
			coupling = 'tri'
				
		self.fix = 'npt temp %f %f %f %s %f %f %f nreset 1000' %(temperature, ramp_to_temp, t_damp, coupling, pressure, pressure, p_damp)
		self.cell_relaxed = True

class Constraint2:
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
		
class Spring(Constraint2):
	def __init__(self, indices, point, spring_constant, R0=0.0):
		Constraint2.__init__(self, indices)
		self.point = point
		self.K = spring_constant
		self.R0 = R0
		
	def get_fix(self, atoms):
		K = atoms.calc.from_ase_units(self.K, 'force')
		x, y, z = atoms.calc.prism.vector_to_lammps(self.point)
		return 'spring tether %f %f %f %f %f' % (K, x, y, z, self.R0)
		
class Constraint:
	def __init__(self):
		self.groups = []
		self.fixes = []
		
	def new_group(self, atoms):
		name = 'group%s%s' %(self.__class__.__name__, abs(hash(tuple(atoms))))
		self.groups.append((name, atoms))
		return name
		
	def get_commands(self, atoms):
		cmds = []
		for name, indices in self.groups:
			indices_str = ' '.join([str(i+1) for i in indices])
			cmds.append('group %s id %s' % (name, indices_str))
		for fix in self.fixes:
			cmds.append(fix)
		return cmds
		
class AddForce(Constraint):
	def __init__(self, target_atoms, total_force):
		Constraint.__init__(self)
		self.target_atoms = target_atoms
		self.total_force = total_force
		
	def get_commands(self, atoms):
		self.fixes = []
		self.groups = []
		
		groupname = self.new_group(self.target_atoms)
		fixname = 'addforce%s' % abs(hash(tuple(self.target_atoms)))
		
		force = self.total_force / len(self.target_atoms)
		fx, fy, fz = atoms.calc.from_ase_units(force, 'force')
		
		self.fixes.append('fix %s %s addforce %s %s %s' % (fixname, groupname, fx, fy, fz))
		return Constraint.get_commands(self, atoms)
		
		