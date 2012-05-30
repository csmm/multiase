from gpaw import GPAW
from ase.io import read 
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.nvtberendsen import NVTBerendsen
from ase import units

restart = 4 #Initial is 0, use integer increment by 1
# Set up the system
if restart == 0:
    from ase.structure import molecule
    atoms = molecule('H2')
    atoms.center(vacuum=2.0)
    atoms.set_pbc(True)
    #atoms = read('relax.traj')
else:
    atoms = read('berendsen_{0}.traj'.format(restart-1))
txt = 'berendsen_{0}.txt'.format(restart)
trajectory = 'berendsen_{0}.traj'.format(restart)
logfile = 'berendsen_{0}.log'.format(restart)

# Describe the interatomic interactions
calc = GPAW(xc='PBE',h=0.3, txt= txt, nbands=-2)
atoms.set_calculator(calc)

# Set the momenta corresponding to T=400K
#MaxwellBoltzmannDistribution(atoms, 0.1*units.kB)

# We want to run MD with constant temperature using Langevin

dyn = NVTBerendsen(atoms, 0.1 * units.fs, 300, taut=100*units.fs,
                   logfile=logfile,
                   trajectory=trajectory)  # 5 fs time step.

dyn.run(20)

