# Copyright (C) 2012 Aalto University
# Author: Lauri Leukkunen <lauri.leukkunen@aalto.fi>

import cPickle as pickle
import numpy as np
from ase.structure import molecule
from ase import Atoms
from ase.constraints import FixBondLengths

methane_count = 500
cell = (100,100,100)
atom_zone = (60,60,60)
atoms = molecule("CH4")
pos_offset = np.array(atom_zone) * np.random.rand(methane_count, 3)
for i in range(methane_count - 1): # the first one is already done
    m = molecule("CH4")
    pos = m.get_positions()
    pos += pos_offset[i]
    m.set_positions(pos)
    atoms.extend(m)

"""
constraints = [0] * methane_count * 4
constraints[0] = [0,1]
constraints[1] = [0,2]
constraints[2] = [0,3]
constraints[3] = [0,4]
for i in range(1,methane_count): # the first one is already done
    m = molecule("CH4")
    pos = m.get_positions()
    pos += pos_offset[i]
    m.set_positions(pos)
    atoms.extend(m)
    constraints[i*4] = [i*5, i*5 + 1]
    constraints[i*4 + 1] = [i*5, i*5 + 2]
    constraints[i*4 + 2] = [i*5, i*5 + 3]
    constraints[i*4 + 3] = [i*5, i*5 + 4]

c = FixBondLengths(constraints)

atoms.center()
atoms.set_constraint(c) # important to set only after setting pos
"""
atoms.center()
fd = open("ch4_gas.pkl", "w+")
pickle.dump(methane_count, fd, -1)
pickle.dump(cell, fd, -1)
pickle.dump(atoms, fd, -1)
fd.close()
