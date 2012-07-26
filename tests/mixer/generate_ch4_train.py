# Copyright (C) 2012 Aalto University
# Author: Lauri Leukkunen <lauri.leukkunen@aalto.fi>

import cPickle as pickle
import numpy as np
from ase.structure import molecule
from ase import Atoms
from ase.constraints import FixBondLengths
from ase.io.trajectory import PickleTrajectory
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase import units

"""
This script generates a cube evenly, but randomly filled out with
CH4 molecules, and two "trains" of CH4 molecules that are set in motion.
"""


atoms = Atoms()
# First do the random molecules
methane_cube_pos = 100.0 * np.random.rand(100, 3) - 50.0

for i in range(len(methane_cube_pos)):
    m = molecule("CH4")
    m.set_positions(m.get_positions() + methane_cube_pos[i])
    atoms.extend(m)

MaxwellBoltzmannDistribution(atoms, 100*units.kB)

# random completed, now the two trains

v = 100.0 / (0.1 * units.fs * 1000) # 1000 simulation steps to cover
                                    # from edge to edge
spread = np.linspace(-50.0, 50.0, 10)
v_array = np.array([v] * (len(spread)))

t1_x = spread
t1_y = np.zeros_like(spread)
t1_z = np.zeros_like(spread)
v1_x = v_array
v1_y = np.zeros_like(v_array)
v1_z = np.zeros_like(v_array)

t2_x = np.zeros_like(spread)
t2_y = spread
t2_z = np.zeros_like(spread)
v2_x = np.zeros_like(v_array)
v2_y = v_array
v2_z = np.zeros_like(v_array)

train1_pos = np.column_stack((t1_x, t1_y, t1_z))
train2_pos = np.column_stack((t2_x, t2_y, t2_z))

train1_vel = np.column_stack((v1_x, v1_y, v1_z))
train2_vel = np.column_stack((v2_x, v2_y, v2_z))

train_pos = np.copy(train1_pos)
train_vel = np.copy(train1_vel)

# second train commented out to avoid collisions that screw GPAW
# calculation
#train_pos = np.append(train_pos, train2_pos, axis=0)
#train_vel = np.append(train_vel, train2_vel, axis=0)

methane_trains = Atoms()

for i in range(len(train_pos)):
    m = molecule("CH4")
    m.set_positions(m.get_positions() + train_pos[i])
    m.set_velocities(np.zeros((5, 3)) + train_vel[i])
    methane_trains.extend(m)

# trains are now ready, append to atoms
#print(methane_trains.get_positions())
atoms.extend(methane_trains)

# now write to .traj file
pt = PickleTrajectory(filename="ch4_train.traj", mode="w", backup=False)
pt.write(atoms)
pt.close()
