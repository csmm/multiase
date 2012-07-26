# Copyright (C) 2012 Aalto University
# Author: Lauri Leukkunen <lauri.leukkunen@aalto.fi>

from copy import deepcopy
import numpy as np

from ase import Atoms
from ase.calculators.interface import Calculator

from csmmcalc.lammps.reaxff import ReaxFF

class Calculation(object):
    """
    Base class for calculations that Mixer can mix.
    You don't want to use this directly, instead look into
    ForceCalculation and EnergyCalculation
    """

    def __init__(self, name=None, selector=None):
        self._name = name
        self._calculator = None
        self._selector = selector
        self.cell = None
        self.coeff = 1.0
 
    def _set_calculator(self, value):
        # check that value at least looks like a calculator
        if not (hasattr(value, "calculation_required") and
                hasattr(value, "get_forces") and
                hasattr(value, "get_potential_energy")):
            raise ValueError("Invalid calculator object given")
        self._calculator = value

    def _get_calculator(self):
        return self._calculator

    def _del_calculator(self):
        del self._calculator

    def _set_selector(self, value):
        self._selector = value
    
    def _get_selector(self):
        return self._selector

    def _del_selector(self):
        del self._selector


    calculator = property(_get_calculator, _set_calculator, 
            _del_calculator, "Calculator property")
    
    selector = property(_get_selector, _set_selector,
                        _del_selector, "Selector property")

    def calculation_required(self, atoms, quantities):
        subset, subset_map, weights = self.get_subset(atoms)
        if not subset:
            return False
        return self.calculator.calculation_required(subset, quantities)

    def get_subset(self, atoms):
        """
        Return the subset of atoms that match the atom_ids of this
        Calculation object. The resulting Atoms object will have
        its atoms centered in the cell defined for this calculation
        and its calculator is also set accordingly.

        @type atoms:           Atoms
        @param atoms:          the normal ASE Atoms object

        @rtype:                (Atoms, dict, list)
        @return:               returns a subset of atoms
                               created with filter_atoms()
                               using self.atom_ids.
                               map[i] gives the original index
                               of Atoms[i] in the atoms list.
        """

        if not self._selector:
            return (atoms, None, {})

        subset, subset_map, weights = self.selector.select_atoms(atoms)
        if not subset:
            return None, None, None
        subset.set_cell(self.cell)
        subset.center()
        subset.set_calculator(self.calculator)
        return subset, subset_map, weights

class ForceCalculation(Calculation):
    """
    ForceCalculation is used by Mixer to drive sub-calculators
    to produce forces for subsets of the atoms. Forces are
    summed using the weights, which can have any float values,
    also negative. If weight is not set for some atom_ids, their
    weight is taken as ForceCalculation.default_weight, which by
    default is 0.0. All atoms listed in atom_ids will take part
    in the calculation irrespective of their associated weight.

    Setting default_weight to 1.0 may be handy to avoid having
    to explicitly declare it for large numbers of atoms in
    case the calculation is the "outer" calculation in a
    multi-scale setup.
    """
    
    def __init__(self, name=None, selector=None, debug=0, debug_file=None):
        super(ForceCalculation, self).__init__(name, selector)
        self.default_weight = 0.0
        self._debug = debug
        self._debug_file = debug_file
        self._debug_counter = 0
        if self._debug > 1 and self._debug_file == None:
            f, fname = tempfile.mkstemp(suffix=".log",
                    prefix="force_calculation_%s-" % self._name,
                    dir=".")
            print("writing debug output to: %s" % fname)
            self._debug_file = io.open(f, mode="w+b", buffering=0)


    def get_forces(self, atoms):
        subset, subset_map, weights = self.get_subset(atoms)
        res = np.zeros((len(atoms), 3))
        if not subset:
            return res
        forces = self.calculator.get_forces(subset)

        if not subset_map:
            return forces
        
        output_forces = self.coeff * weights * forces
        
        if self._debug > 1:
            self._debug_counter += 1
            atom_ids = Mixer.get_atom_ids(subset)
            for i in range(len(forces)):
                self._debug_file.write(
                    "%s,%i,%f,%s,%s,%s" % (
                        self._name,
                        atom_ids[i],
                        self._debug_counter,
                        self.coeff, weights[i], forces[i],
                        output_forces[i]))

        
        # now map them to the original atom sequence
        for i in range(len(subset)):
            res[subset_map[i]] = output_forces[i]
        
        return res


class EnergyCalculation(Calculation):
    """
    EnergyCalculation is used by Mixer to drive
    sub-calculators to produce energies for subsets
    of the atoms.
    """
    def __init__(self, name=None, selector=None):
        super(EnergyCalculation, self).__init__(name, selector)

    def get_energy(self, atoms, force_consistent=False):
        """
        Calculates the energy for the subset of atoms
        defined for this EnergyCalculation

        @type atoms:               Atoms
        @param atoms:               ASE Atoms object
        @type force_consistent:     Boolean
        @param force_consistent:    (ignored, False assumed always)
        """
        subset, subset_map, weights = self.get_subset(atoms)
        if not subset:
            return 0.0
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
    def __init__(self, forces, energies, debug=0):
        """
        @type forces:           [ForceCalculation]
        @param forces:          list of Mixer.ForceCalculation to 
                                drive force calculations
        @type energies:         [EnergyCalculation]
        @param energies:        list of Mixer.EnergyCalculation to drive 
                                desired number of energy region
                                calculations
        """

        # check that forces and energies are not obviously broken
        
        for c in forces + energies:
            if (c.calculator == None or
                c.cell == None):
                raise ValueError("Invalid/Uninitialized Calculation" +
                    " object passed")

        self._forces = forces
        self._energies = energies
        self._debug = debug

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
        return list(atoms.get_array("csmmcalc.mixer.atom_ids"))

