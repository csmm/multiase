from ase.optimize.optimize import Dynamics

class LAMMPSOptimizer(Dynamics):
	""" Geometry optimizer for LAMMPS. works only with LAMMPS calculators """
	def __init__(self, atoms, restart=None, logfile=None, trajectory=None):
		Dynamics.__init__(self, atoms, logfile, trajectory)
		
	def run(self, fmax=0.05, steps=1e8):
		positions = self.atoms.calc.minimize(self.atoms, ftol=fmax, maxeval=steps)
		self.atoms.positions = positions

		
class LAMMPSMolecularDynamics(Dynamics):
	""" Base class for molecular dynamics with LAMMPS. Requires a LAMMPS calculator. """
	def __init__(self, atoms, timestep, integrator='verlet',
				trajectory=None, logfile=None, restart=None, loginterval=1):
		Dynamics.__init__(self, atoms, logfile, trajectory)
		
		self.dt = timestep
		
		if integrator == 'verlet':
			self.run_style = 'velet'
		else:
			raise RuntimeError('Unknown integrator: %s' % thermostat)
		
		self.fix = None
	
		
	def run(self, steps=50):
		positions, momenta = self.atoms.calc.molecular_dynamics(
									self.atoms, self.dt, steps, 'all '+self.fix)
		self.atoms.positions = positions
		self.atoms.set_momenta(momenta)
		

class LAMMPS_NVE(LAMMPSMolecularDynamics):
	""" Microcanonical ensemble """
	def __init__(self, atoms, timestep, integrator='verlet',
				trajectory=None, logfile=None, loginterval=1):
		LAMMPSMolecularDynamics.__init__(self, atoms, timestep, integrator, trajectory, logfile, loginterval)
		
		self.fix = 'nve'
		
class LAMMPS_NVT(LAMMPSMolecularDynamics):
	""" Constant temperature calculations with Nose-Hoover or Berendsen """
	def __init__(self, atoms, timestep, temperature, t_damp=100, thermostat='Nose-Hoover',
				integrator='verlet', trajectory=None, logfile=None, loginterval=1):
		LAMMPSMolecularDynamics.__init__(self, atoms, timestep, integrator, trajectory, logfile, loginterval)
		
		if thermostat == 'Nose-Hoover':
			cmd = 'nvt temp'
		elif thermostat == 'Berendsen':
			cmd = 'temp/berendsen'
		else:
			raise RuntimeError('Unknown thermostat: %s' % thermostat)
		
		self.fix = '%s %f %f %f' %(cmd, temperature, temperature, t_damp)