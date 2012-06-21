# Copyright (C) 2012 Aalto University
# Author: lauri.leukkunen@aalto.fi

from copy import deepcopy
import numpy as np

from ase import Atoms
from ase.calculators.interface import Calculator

from gpaw import GPAW

from csmmcalc.lammps.reaxff import ReaxFF

class Mixer(Calculator):
    def __init__(self, calcs):
        """
        @param calcs    list of initialized ASE calculator objects
        """
        self._atoms = [0] * len(calcs)
        self._calcs = calcs


    def calculation_required(self, atoms, quantities):
        for c in self._calcs:
            if c.calculation_required(atoms, quantities):
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
        return None

    def _get_atoms(self, atoms):
        """
        @param atoms    - the normal ASE Atoms object

        @return (atoms[], weights[], maps[])
            - returns a tuple that contains three lists that
              have the atoms, mixing weights and atom mappings
              for each calculator
        """
        # create per calc atom lists
        atom_lists = [[] for x in range(len(self._calcs))]
        weight_lists = [[] for x in range(len(self._calcs))]
        my_atoms = [0] * len(self._calcs)
        my_weights = [0] * len(self._calcs)
        my_maps = [[] for x in range(len(self._calcs))]
        
        for i in range(len(atoms)):
            t = atoms[i].tag
            for j in range(len(self._calcs)):
                weight = Mixer.weight(t, j)
                if weight > 0.0:
                    atom_lists[j].append(atoms[i:i+1])
                    weight_lists[j].append(weight)
                    my_maps[j].append(i)

        for i in range(len(self._calcs)):
            my_atoms[i] = Atoms(atom_lists[i][0])
            for a in atom_lists[i][1:]:
                my_atoms[i].extend(a)
            my_weights[i] = weight_lists[i]
        
        return (my_atoms, my_weights, my_maps)
            
    @staticmethod
    def tag(calc, weight):
        """
        @param calc     calculator number (index starts from 0)
        @param weight   float in range 0.0 - 1.0, This is mapped to a 16 bit
                        int.

        @return int
        """
        if calc < 0 or (weight < 0.0 or weight > 1.0):
            raise ValueError("Unable to create calculator weight")

        return int(0xFFFF & int(weight * 0xFFFF)) << calc * 16

    @staticmethod
    def weight(tag, calc):
        """
        @param tag              ASE tag encoded using tag_weight()
        @param calc             which calc weight to extract

        @return float           returns float in range 0.0 - 1.0
        """
        return float(float(0xFFFF & (tag >> calc * 16)) / float(0xFFFF))

