from gpaw import GPAW
from ase.io import read 
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units

restart = 6 #Initial is 0, use integer increment by 1
# Set up the system
if restart == 0:
    from ase.structure import molecule
    atoms = molecule('H2')
    atoms.center(vacuum=2.0)
    atoms.set_pbc(True)
    #atoms = read('relax.traj')
else:
    atoms = read('verlet_{0}.traj'.format(restart-1))
txt = 'verlet_{0}.txt'.format(restart)
trajectory = 'verlet_{0}.traj'.format(restart)
logfile = 'verlet_{0}.log'.format(restart)

# Describe the interatomic interactions
calc = GPAW(xc='PBE',h=0.3, txt= txt,
            nbands=-2)
atoms.set_calculator(calc)

# Set the momenta corresponding to T=400K
if restart == 0:
    MaxwellBoltzmannDistribution(atoms, 400*units.kB)

# We want to run MD with constant energy using the VelocityVerlet algorithm.
dyn = VelocityVerlet(atoms, 1*units.fs, logfile=logfile,
                     trajectory=trajectory)  # 5 fs time step.

dyn.run(5)
    
