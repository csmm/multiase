# Copyright (C) 2012 Aalto University
# Author: Lauri Leukkunen <lauri.leukkunen@aalto.fi>

import cPickle as pickle
import numpy as np
from ase.structure import molecule
from ase import Atoms
from ase.constraints import FixBondLengths

methane_count = 10
edge = 100.0
cell = (edge,edge,edge)
atom_zone = cell
atoms = molecule("CH4")
atoms.set_positions(atoms.get_positions() + edge/2.0)
pos_offset = np.array(atom_zone) * np.random.rand(methane_count, 3)
for i in range(methane_count - 1): # the first one is already done
    m = molecule("CH4")
    pos = m.get_positions()
    pos += pos_offset[i]
    m.set_positions(pos)
    atoms.extend(m)

atoms.set_positions(atoms.get_positions() - edge/2.0)
fd = open("ch4_gas.pkl", "w+")
pickle.dump(methane_count, fd, -1)
pickle.dump(cell, fd, -1)
pickle.dump(atoms, fd, -1)
fd.close()
