# Copyright (C) 2012 Aalto University
# Author: lauri.leukkunen@aalto.fi

from copy import deepcopy

from ase.calculators.interface import Calculator

from csmmcalc.lammps.reaxff import ReaxFF

class Mixer(Calculator):
    def calculation_required(self, atoms, quantities):
        return False

    def get_forces(self, atoms):
        return None

    def get_potential_energy(self, atoms=None, force_consistent=False):
        return None

    def get_stress(self, atoms):
        return None

    def set_atoms(self, atoms):
        self._atoms = copy(atoms)

