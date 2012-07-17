# Copyright (C) 2012 Aalto University
# Author: Lauri Leukkunen <lauri.leukkunen@aalto.fi>

import cPickle as pickle
import numpy as np
from ase.structure import molecule
from ase import Atoms
from ase.constraints import FixBondLengths

methane_count = 20
edge = 100.0
cell = (edge,edge,edge)

atoms = Atoms()
pos_offset = None
if 0:
    pos_offset = np.array(atom_zone) * np.random.rand(methane_count, 3)
else:
    a = np.zeros((methane_count, 3))
    step = (edge-20)/(methane_count-1)
    for i in range(methane_count):
        a[i][0] = -(edge - 20)/2 + i * step
    pos_offset = a

for i in range(methane_count):
    m = molecule("CH4")
    pos = m.get_positions()
    pos += pos_offset[i]
    m.set_positions(pos)
    atoms.extend(m)

print(atoms.get_positions())
fd = open("ch4_gas.pkl", "w+")
pickle.dump(methane_count, fd, -1)
pickle.dump(cell, fd, -1)
pickle.dump(atoms, fd, -1)
fd.close()
