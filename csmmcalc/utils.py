import os
from os import path
from ase import units
from ase.md.md import MolecularDynamics
import numpy as np


csmm_config = {
        "CSMM_INSTALL_DIR":path.split(path.dirname(path.realpath(__file__)))[0],
        }

def get_config(key):
    value = os.getenv(key)
    if not value:
        value = csmm_config[key]
    return value

def get_datafile(filename):
    datadir = os.path.join(get_config("CSMM_INSTALL_DIR"), "data")
    return os.path.join(datadir, os.path.basename(filename))


class DynTesting(MolecularDynamics):
    def __init__(self, atoms, timestep=0.1*units.fs,
                 logfile=None, trajectory=None,
                 loginterval=1, offset=(0.1, 0, 0)):

        MolecularDynamics.__init__(self, atoms, timestep,
                                   logfile, trajectory,
                                   loginterval)
        self._offset = np.array(offset)

    def step(self, f):
        print("stepping")
        f = self.atoms.get_forces()
        self.atoms.set_positions(self.atoms.get_positions() + self._offset)
        return f

