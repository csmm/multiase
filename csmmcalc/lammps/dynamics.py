from ase.optimize.optimize import Dynamics
from ase.io.trajectory import PickleTrajectory
from ase.md.logger import MDLogger
from ase import units

class LAMMPSOptimizer(Dynamics):
	""" Geometry optimizer for LAMMPS. works only with LAMMPS calculators """
	def __init__(self, atoms, restart=None, logfile=None, trajectory=None):
		Dynamics.__init__(self, atoms, logfile, trajectory)
		
	def run(self, fmax=0.001, steps=1e8, min_style=None):
		self.atoms.calc.minimize(self.atoms, ftol=fmax, maxeval=steps, min_style=min_style)

		
class LAMMPSMolecularDynamics(Dynamics):
	""" Base class for molecular dynamics with LAMMPS. Requires a LAMMPS calculator. """
	def __init__(self, atoms, timestep, integrator='verlet',
				trajectory=None, traj_interval=1000, logfile=None, restart=None, loginterval=100):
		Dynamics.__init__(self, atoms, None, None)

		self.dt = timestep
		
		if integrator == 'verlet':
			self.run_style = 'velet'
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
		
	def run(self, steps=50):
		self.nsteps = 0
		fix = 'all '+self.fix
		calculation = self.atoms.calc.molecular_dynamics(self.atoms, self.dt, fix)
		calculation.next()
		for target_step in xrange(0, steps+1):
			for function, interval, args, kwargs in self.observers:
				if target_step % interval == 0:
					if target_step > self.nsteps:
						# Calculate up to nsteps
						calculation.send(target_step - self.nsteps)
						self.nsteps = target_step
					function(*args, **kwargs)
					
		if target_step > self.nsteps:
			calculation.send(target_step - self.nsteps)
		try:
			calculation.send(None)
		except StopIteration:
			pass
	
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
				thermostat='Nose-Hoover', **kwargs):
		LAMMPSMolecularDynamics.__init__(self, atoms, timestep, **kwargs)
		
		if thermostat == 'Nose-Hoover':
			cmd = 'nvt temp'
		elif thermostat == 'Berendsen':
			cmd = 'temp/berendsen'
		else:
			raise RuntimeError('Unknown thermostat: %s' % thermostat)
		
		t_damp = atoms.calc.from_ase_units(t_damp, 'time')
		
		self.fix = '%s %f %f %f' %(cmd, temperature, temperature, t_damp)