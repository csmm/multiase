# Copyright (C) 2012 Aalto University
# Author: lauri.leukkunen@aalto.fi

from copy import deepcopy
import numpy as np

from ase import Atoms
from ase.calculators.interface import Calculator

from gpaw import GPAW

from csmmcalc.lammps.reaxff import ReaxFF


class Calculation(object):
    calculator = None
    cell = None
    atom_tags = ()

    def calculation_required(self, atoms, quantities):
        subset, subset_map = self.get_subset(atoms)
        return self.calculator.calculation_required(subset, quantities)

    def get_subset(self, atoms):
        subset, subset_map = self.filter_atoms(atoms, self.atom_tags)
        subset.set_cell(self.cell)
        subset.center()
        subset.set_calculator(self.calculator)
        return subset, subset_map

    def filter_atoms(self, atoms, tags):
        """
        @param atoms    the normal ASE Atoms object
        @param tags     list of tags to match

        @return Atoms, map
                        map[i] gives the original index
                        of Atoms[i]
        """
        subset = None
        subset_map = {}
        k = 0
        for i in range(len(atoms)):
            if atoms[i].tag not in tags:
                continue

            if not subset:
                subset = Atoms(atoms[i:i+1])
            else:
                subset.extend(atoms[i:i+1])

            subset_map[k] = i
            k += 1

        return subset, subset_map


class ForceCalculation(Calculation):
    weights = ()
    
    def get_forces(self, atoms):
        subset, subset_map = self.get_subset(atoms)
        
        forces = self.calculator.get_forces(subset)

        weight_array = np.ones_like(forces)
        for i in range(len(subset)):
            weight_array[i] = weight_array[i] * self.weights[self.atom_tags[i]]
        
        forces *= weight_array
        
        # now map them to the original atom sequence
        res = np.zeros((len(atoms), 3))
        for i in range(len(subset)):
            res[subset_map[i]] = forces[i]
        
        return res


class EnergyCalculation(Calculation):
    coeff = 1

    def get_energy(self, atoms, force_consistent=False):
        subset, subset_map = self.get_subset(atoms)
        energy = self.calculator.get_potential_energy(subset)
        return self.coeff * energy


class Mixer(Calculator):
    def __init__(self, forces=None, energies=None):
        """
        @param forces           list of Mixer.ForceCalculation to drive force
                                calculations
        @param energies         list of Mixer.EnergyCalculation to drive desired
                                number of energy region calculations
        """
        self._forces = forces
        self._energies = energies

    def calculation_required(self, atoms, quantities):
        for c in self._forces + self._energies:
            if c.calculation_required(atoms, quantities):
                return True
        return False
        

    def get_forces(self, atoms):
        """
        @param atoms        the normal ASE Atoms object
        @return np.array()  the expected np.array() containing
                            the forces
        """
        forces = np.zeros((len(atoms), 3))
        for fc in self._forces:
            forces += fc.get_forces(atoms)
        return forces

    def get_potential_energy(self, atoms=None, force_consistent=False):
        """
        @param atoms        the normal ASE Atoms object
        @return float       the energy of the system
        """
        energy = 0.0
        for ec in self._energies:
            energy += ec.get_energy(atoms, force_consistent)
        return energy

