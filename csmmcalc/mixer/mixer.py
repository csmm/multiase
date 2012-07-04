# Copyright (C) 2012 Aalto University
# Author: lauri.leukkunen@aalto.fi

from copy import deepcopy
import numpy as np

from ase import Atoms
from ase.calculators.interface import Calculator

from gpaw import GPAW

from csmmcalc.lammps.reaxff import ReaxFF


class Calculation(object):
    """
    Base class for calculations that Mixer can mix.
    You don't want to use this directly, instead look into
    ForceCalculation and EnergyCalculation
    """
    calculator = None
    cell = None
    atom_ids = ()

    def calculation_required(self, atoms, quantities):
        subset, subset_map = self.get_subset(atoms)
        return self.calculator.calculation_required(subset, quantities)

    def get_subset(self, atoms):
        """
        Return the subset of atoms that match the atom_ids of this
        Calculation object.

        @type atoms:           Atoms
        @param atoms:          the normal ASE Atoms object

        @rtype:                (Atoms, dict)
        @return:               returns a subset of atoms
                               created with filter_atoms()
                               using self.atom_ids.
                               map[i] gives the original index
                               of Atoms[i] in the atoms list.
        """
        subset, subset_map = self.filter_atoms(atoms, self.atom_ids)
        subset.set_cell(self.cell)
        subset.center()
        subset.set_calculator(self.calculator)
        return subset, subset_map

    def filter_atoms(self, atoms, atom_ids):
        """
        Return an Atoms object containing only those from atoms than
        match the given atom_ids list.

        @param atoms:       the normal ASE Atoms object
        @param atom_ids:    list of atom ids to match

        @rtype:             (Atoms, dict)
        @return:            map[i] gives the original index
                            of Atoms[i]
        """
        subset = None
        subset_map = {}
        k = 0
        ids = Mixer.get_atom_ids(atoms)
        for i in range(len(atoms)):
            if ids[i] not in atom_ids:
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
            weight_array[i] = weight_array[i] * self.weights[self.atom_ids[i]]
        
        forces *= weight_array
        
        # now map them to the original atom sequence
        res = np.zeros((len(atoms), 3))
        for i in range(len(subset)):
            res[subset_map[i]] = forces[i]
        
        return res


class EnergyCalculation(Calculation):
    """
    EnergyCalculation is used by Mixer to drive
    sub-calculators to produce energies for subsets
    of the atoms.
    """
    coeff = 1

    def get_energy(self, atoms, force_consistent=False):
        """
        Calculates the energy for the subset of atoms
        defined for this EnergyCalculation

        @type atoms:               Atoms
        @param atoms:               ASE Atoms object
        @type force_consistent:     Boolean
        @param force_consistent:    (ignored, False assumed always)
        """
        subset, subset_map = self.get_subset(atoms)
        energy = self.calculator.get_potential_energy(subset)
        return self.coeff * energy


class Mixer(Calculator):
    """
    Mixer allows combining any number of real ASE calculators.
    Typically this is used in a multi-scale calculation where
    different levels of theory are used for different parts of
    the system. One example would be combining GPAW with ReaxFF.

    Please look at tests/mixer/h2_mixer.py for an example of how
    to use this calculator.
    """
    def __init__(self, forces, energies):
        """
        @type forces:           [ForceCalculation]
        @param forces:          list of Mixer.ForceCalculation to drive force
                                calculations
        @type energies:         [EnergyCalculation]
        @param energies:        list of Mixer.EnergyCalculation to drive desired
                                number of energy region calculations
        """

        # check that forces and energies are not obviously broken
        
        for c in forces + energies:
            if (c.calculator == None or
                c.atom_ids == None or
                c.cell == None):
                raise ValueError("Invalid/Uninitialized Calculation object passed")

        self._forces = forces
        self._energies = energies

    def calculation_required(self, atoms, quantities):
        for c in self._forces + self._energies:
            if c.calculation_required(atoms, quantities):
                return True
        return False
        

    def get_forces(self, atoms):
        """
        @type atoms:        Atoms
        @param atoms:       the normal ASE Atoms object

        @rtype:             np.array
        @return:            the expected np.array containing
                            the forces
        """
        forces = np.zeros((len(atoms), 3))
        for fc in self._forces:
            forces += fc.get_forces(atoms)
        return forces

    def get_potential_energy(self, atoms=None, force_consistent=False):
        """
        @type atoms:        Atoms
        @param atoms:       the normal ASE Atoms object

        @rtype:             float
        @return:            the energy of the system
        """
        energy = 0.0
        for ec in self._energies:
            energy += ec.get_energy(atoms, force_consistent)
        return energy
    
    def set_atoms(self, atoms):
        self.atoms = atoms.copy() # this is done to comply with
                                  # the ASE interface.py

        # create id array if it is not yet set
        if not atoms.has("csmmcalc.mixer.atom_ids"):
            Mixer.set_atom_ids(atoms)

    @staticmethod
    def set_atom_ids(atoms):
        ids = np.array(range(len(atoms)))
        atoms.new_array("csmmcalc.mixer.atom_ids", ids)

    @staticmethod
    def get_atom_ids(atoms):
        return atoms.get_array("csmmcalc.mixer.atom_ids")
