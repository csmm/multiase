# Copyright (C) 2012 Aalto University
# Author: lauri.leukkunen@aalto.fi

from copy import deepcopy
import numpy as np

from ase import Atoms
from ase.calculators.interface import Calculator

from gpaw import GPAW

from csmmcalc.lammps.reaxff import ReaxFF

class Mixer(Calculator):
    def __init__(self, calcs, cells=None):
        """
        @param calcs        list of initialized ASE calculator objects
        @param cells        list of unit cell matrices to use, one per
                            calculator
        """
        self._atoms = [0] * len(calcs)
        self._calcs = calcs
        self._cells = cells


    def calculation_required(self, atoms, quantities):
        a,w,m = self._get_atoms(atoms)
        for i in range(len(self._calcs)):
            if self._calcs[i].calculation_required(a[i], quantities):
                return True
        return False
        

    def get_forces(self, atoms):
        """
        @param atoms    - the normal ASE Atoms object
        @return array() - the expected np.array() containing
                          the forces
        """
        forces = np.zeros((len(atoms), 3))
        a_list, w_list, a_maps = self._get_atoms(atoms)
        for i in range(len(self._calcs)):
            f = self._calcs[i].get_forces(a_list[i])
            wa = np.ones_like(f)
            for j in range(len(w_list[i])):
                wa[j] = wa[j] * w_list[i][j]
            f = f * wa
            for k in range(len(a_list[i])):
                forces[a_maps[i][k]] += f[k]

        return forces

    def get_potential_energy(self, atoms=None, force_consistent=False):
        """
        TODO: What should be done here? The mixed hamiltonian model
              doesn't really work for potential energy calculations
              if we don't get contributions of each atom separately.

              We *could* ask the sub-calcs to calculate energy for
              configuration where each atom with 0<w<1.0 is left out
              one by one and thus get their individual contribution
              and mix using the same weights as with forces?
        """
        return 0.0

    def _get_atoms(self, atoms):
        """
        @param atoms    - the normal ASE Atoms object

        @return (atoms[], weights[], maps[])
            - returns a tuple that contains three lists that
              have the atoms, mixing weights and atom mappings
              for each calculator
        """
        # create per calc atoms, weight, map lists
        
        wrk_atoms = [0] * len(self._calcs)
        wrk_weights = [[] for x in range(len(self._calcs))]
        wrk_maps = [[] for x in range(len(self._calcs))]
        
        for i in range(len(atoms)):
            t = atoms[i].tag
            for j in range(len(self._calcs)):
                weight = Mixer.weight(j, t)
                if weight > 0.0:
                    if wrk_atoms[j] == 0:
                        wrk_atoms[j] = Atoms(atoms[i:i+1])
                    else:
                        wrk_atoms[j].extend(atoms[i:i+1])
                    wrk_weights[j].append(weight)
                    wrk_maps[j].append(i)

        # Atoms object of each calculator now contain the right atoms.
        # All that's left is to switch to right unit cell and center
        # the atoms so that they hopefully fit within it.

        for i in range(len(self._calcs)):
            wrk_atoms[i].set_cell(self._cells[i])
            wrk_atoms[i].center()

        return (wrk_atoms, wrk_weights, wrk_maps)
            

    BITS_PER_WEIGHT = 16

    @staticmethod
    def tag(calc, weight):
        """
        @param calc     calculator number (index starts from 0)
        @param weight   float in range 0.0 - 1.0, This is mapped to an
                        int of Mixer.BITS_PER_WEIGHT bits (typically
                        32).

        @return int
        """
        if calc < 0 or (weight < 0.0 or weight > 1.0):
            raise ValueError("Unable to create calculator weight")

        return int(0xFFFF & int(weight * 0xFFFF)) << (calc *
            Mixer.BITS_PER_WEIGHT)

    @staticmethod
    def weight(calc, tag):
        """
        @param calc             which calc weight to extract
        @param tag              ASE tag encoded using Mixer.tag()

        @return float           returns float in range 0.0 - 1.0
        """
        return float(float(0xFFFF & (tag >> calc *
            Mixer.BITS_PER_WEIGHT)) / float(0xFFFF))

