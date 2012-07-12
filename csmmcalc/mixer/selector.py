# Copyright (C) 2012 Aalto University
# Author: Lauri Leukkunen <lauri.leukkunen@aalto.fi>


from ase import Atoms
import numpy as np

from mixer import Mixer

class AtomSelector(object):
    """
    Base class for all atom selectors. Purpose of these selectors
    is to provide Atoms objects for ASE calculators to use within 
    the context of Mixer calculations.
    """
    def select_atoms(self, atoms):
        """
        Return an Atoms object containing only those from atoms than
        match the atom_ids list.

        @type atoms:        Atoms
        @param atoms:       the normal ASE Atoms object
        @type atom_ids:     list
        @param atom_ids:    list of atom ids to match

        @rtype:             (Atoms, dict, list)
        @return:            Atoms provides an ASE Atoms object with the selected
                            atoms. Dict contains a map, where map[i] gives the
                            original index of Atoms[i]. List provides force
                            calculation weights for each atom in Atoms.
        """
        raise NotImplementedError("Override this")

class AtomListSelector(AtomSelector):
    """
    Very simple atom selector that selects based on a list of
    atom ids. This can be used if you know the atoms do not move around
    too much and want to fix them permanently to a calculation. Naturally
    this list can be updated during a long multi-step computation, but
    for a selector that provides a better spatial calculation region support,
    see CalcBox.
    """
    def __init__(self, atom_ids=None, weights=None):
        super(AtomListSelector, self).__init__()
        self._atom_ids = atom_ids
        self._weights = weights

    def set_atom_ids(self, atom_ids):
        """
        Set the atom ids list used for selecting.

        @type atom_ids:         list
        @param atom_ids:        list of atom_ids
        """
        self._atom_ids = atom_ids

    def set_weights(self, weights):
        """
        Set the force calculation weights.

        @type weights:          dict
        @param weights:         force calculation weights
        """
        self._weights = weights

    def select_atoms(self, atoms):
        subset = None
        subset_map = {}
        subset_weights = []
        k = 0
        ids = Mixer.get_atom_ids(atoms)
        for i in range(len(atoms)):
            if ids[i] not in self._atom_ids:
                continue

            if not subset:
                subset = Atoms(atoms[i:i+1])
            else:
                subset.extend(atoms[i:i+1])

            subset_map[k] = i
            subset_weights.append(self._weights[ids[i]])
            k += 1
        wa = np.zeros((len(subset_weights), 3))
        for i in range(len(subset_weights)):
            wa[i] += subset_weights[i]
        return subset, subset_map, wa



class CalcRegion(AtomSelector):
    """
    Base class for spatially oriented selectors. These can be used
    to carve up the system space for different calculators while
    allowing atoms to move around freely.
    """
    def __init__(self, pos=(0,0,0), pbc=None):
        """
        @type pos:      tuple(float, float, float)
        @param pos:     Coordinates for the center of the region
        
        """
        super(CalcRegion, self).__init__()
        self._pos = pos
        self._pbc = None

    def atom_inside(self, pos):
        """
        Tests if given coordinates are within the CalcRegion.

        @type pos:          tuple(float, float, float)
        @param pos:         coordinates of the atom

        @rtype:             Boolean
        @return:            True if inside, False if not
        """
        raise NotImplementedError("Override this")

    def get_weight(self, pos):
        """
        Returns a weight for the coordinates.

        @type pos:          tuple(float, float, float)
        @param pos:         coordinates to get a weight for

        @rtype:             float
        @return:            The weight.
        """
        raise NotImplementedError("Override this")

    def select_atoms(self, atoms):
        subset = None
        subset_map = {}
        weights = []
        k = 0
        for i in range(len(atoms)):
            if self.atom_inside(atoms[i].position):
                if not subset:
                    subset = atoms[i:i+1]
                else:
                    subset.extend(atoms[i:i+1])
                weights.append(self.get_weight(atoms[i].position))
            subset_map[k] = i
            k += 1

        wa = np.zeros((len(weights), 3))
        for i in range(len(weights)):
            wa[i] += weights[i]
        
        if self._pbc:
            subset.set_pbc(self._pbc)
        return subset, subset_map, wa


class CalcBox(CalcRegion):
    """
    Atom selector for rectangular volumes with optional transition buffer
    region. The transition region adjusts weights starting from 0.0 on the outer
    surface and reaching 1.0 on the inner surface. The volume contained within
    the inner surface gives weight 1.0.

    The volume can be made periodic by setting the pbc parameter for the
    constructor. Its simply passed on to the generated Atoms object, so it
    should adhere to the ASE documentation.
    """
    def __init__(self, pos=(0,0,0), dim=(1,1,1), inner_dim=None, pbc=None):
        """
        @type pos:          tuple(float, float, float)
        @param pos:         center of the calculation box
        @type dim:          tuple(float, float, float)
        @param dim:         outer dimensions of the box
        @type inner_dim:    tuple(float, float, float)
        @param inner_dim:   dimensions of the full-weight inner box
        @type pbc:          (Boolean, Boolean, Boolean) or Boolean
        @param pbc:         periodic boundary condition flags
        """

        super(CalcBox, self).__init__(pos, pbc)

        self._xmin = pos[0] - dim[0]/2.0
        self._xmax = pos[0] + dim[0]/2.0
        self._ymin = pos[1] - dim[1]/2.0
        self._ymax = pos[1] + dim[1]/2.0
        self._zmin = pos[2] - dim[2]/2.0
        self._zmax = pos[2] + dim[2]/2.0

        self._inner_dim = inner_dim

        if inner_dim:
            self._inner_xmin = pos[0] - inner_dim[0]/2.0
            self._inner_xmax = pos[0] + inner_dim[0]/2.0
            self._inner_ymin = pos[0] - inner_dim[1]/2.0
            self._inner_ymax = pos[0] + inner_dim[1]/2.0
            self._inner_zmin = pos[0] - inner_dim[2]/2.0
            self._inner_zmax = pos[0] + inner_dim[2]/2.0

    def get_weight(self, pos):
        x, y, z = pos
        if (x > self._xmax or x < self._xmin
                or y > self._ymax or y < self._ymin
                or z > self._zmax or z < self._zmin):
            return 0.0

        if not self._inner_dim:
            return 1.0
        
        r = []

        if x > self._xmin and x < self._inner_xmin:
            d = self._inner_xmin - x
            w = self._scale(d, self._xmin, self._inner_xmin)
            r.append((d, w))
        if x < self._xmax and x > self._inner_xmax:
            d = x - self._inner_xmax
            w = self._scale(d, self._inner_xmax, self._xmax)
            r.append((d, w))
        if y > self._ymin and y < self._inner_ymin:
            d = self._inner_ymin - y
            w = self._scale(d, self._ymin, self._inner_ymin)
            r.append((d, w))
        if y < self._ymax and y > self._inner_ymax:
            d = y - self._inner_ymax
            w = self._scale(d, self._inner_ymax, self._ymax)
            r.append((d, w))
        if z > self._zmin and z < self._inner_zmin:
            d = self._inner_zmin - z
            w = self._scale(d, self._zmin, self._inner_zmin)
            r.append((d, w))
        if z < self._zmax and z > self._inner_zmax:
            d = z - self._inner_zmax
            w = self._scale(d, self._inner_zmax, self._zmax)
            r.append((d, w))

        if len(r) == 0:
            return 1.0

        weight = 0.0

        for d, w in r:
            weight = max(weight, w)

        return weight

    @staticmethod
    def _scale(value, left, right):
        return (value - left)/(right - left)


    def atom_inside(self, pos):
        if (pos[0] >= self._xmin and
            pos[0] <= self._xmax and
            pos[1] >= self._ymin and
            pos[1] <= self._ymax and
            pos[2] >= self._zmin and
            pos[2] <= self._zmax):
            return True
        else:
            return False


