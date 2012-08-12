# Copyright (C) 2012 Aalto University
# Author: Lauri Leukkunen <lauri.leukkunen@aalto.fi>

import tempfile
import os
import io

from ase import Atoms
from ase.io.trajectory import PickleTrajectory
import ase.io.pdb

import numpy as np

from mixer import Mixer
from octree import OctreeNode

class AtomSelector(object):
    """
    Base class for all atom selectors. Purpose of these selectors
    is to provide Atoms objects for ASE calculators to use within 
    the context of Mixer calculations.
    """
    def select_atoms(self, atoms):
        """
        Return an Atoms object containing only those from atoms that
        match the selection criteria of the AtomSelector object.

        @type atoms:        Atoms
        @param atoms:       the normal ASE Atoms object

        @rtype:             (Atoms, dict, list)
        @return:            Atoms provides an ASE Atoms object with
                            the selected atoms. Dict contains a map,
                            where map[i] gives the original index of
                            Atoms[i]. List provides force calculation
                            weights for each atom in Atoms.
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
    def __init__(self, atom_ids, weights):
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
    def __init__(self, name, pos=(0,0,0), cutoff=2.0, pbc=None, debug=0):
        """
        @type pos:      tuple(float, float, float)
        @param pos:     Coordinates for the center of the region
        
        """
        super(CalcRegion, self).__init__()
        self._pos = pos
        self._cutoff = cutoff
        self._pbc = None
        self._debug = debug
        self._debug_file = None
        self._name = name
        if self._debug > 1:
            fname = "calcregion_debug_%s.pdb" % self._name
            self._debug_file = open(fname, "w+b", buffering=0)
            print("writing debug output to: %s" % fname)

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

    def get_bounding_box(self):
        """
        Returns a rectangular bounding box for the CalcRegion.
        
        @rtype:             tuple(position,dimensions)
        @return:            Position and dimensions are both
                            tuples of three floats. Position
                            should be the center of the CalcRegion
                            and dimensions the total lengths of each side.
        """
        raise NotImplementedError("Override this")

    def length2(self, a, b):
        """
        Utility function for getting square of distance between
        two points.

        @type a:            tuple(float, float, float)
        @param a:           first coordinate
        @type b:            tuple(float, float, float)
        @param b:           second coordinate

        @rtype:             float
        @return:            distance between points a and b
        """
        return (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2

    def select_atoms(self, atoms):
        """
        Return an Atoms object containing only those from atoms that
        are within the CalcRegion and that don't have an interaction-
        chain to an atom outside the CalcRegion. This interaction-chain
        is a chain of atoms within cutoff distance from each other.

        @type atoms:        Atoms
        @param atoms:       the normal ASE Atoms object

        @rtype:             (Atoms, dict, list)
        @return:            Atoms provides an ASE Atoms object with the
                            selected atoms. Dict contains a map, where map[i]
                            gives the original index of Atoms[i]. List provides
                            force calculation weights for each atom in Atoms.
        """
        subset = None
        subset_map = {}
        weights = []
        cutoff2 = self._cutoff**2 # use square for comparison with length2()
        
        p, d = self.get_bounding_box() # get the rectangular bounding box

        # create Octree with atoms inside or within cutoff distance from
        # CalcRegion
        root = OctreeNode(pos=p,
                dim=(d[0] + 2*self._cutoff,
                     d[1] + 2*self._cutoff,
                     d[2] + 2*self._cutoff),
                res=max(self._cutoff/2.0, 0.1))

        atom_ids = Mixer.get_atom_ids(atoms) # extract atomids from atoms
        rev_map = {} # reverse map used to get from atom_ids to atoms[] indexes
        
        for i in range(len(atoms)):
            root.add_object(atom_ids[i], atoms[i].position)
            rev_map[atom_ids[i]] = i
        
        black_list = {} # blacklist of atoms with interaction-chain to outside
        for k in root.get_objects():
            black_list[k] = False

        wrk = [] # stack for holding blacklisted atoms whose neighbours
                 # need to be checked 

        # populate wrk with those atoms outside CalcRegion that are still
        # within cutoff distance from it
        for atom_id in root.get_objects():
            if self.atom_inside(atoms[rev_map[atom_id]].position):
                continue
            wrk.append(atom_id)
            black_list[atom_id] = True

        # loop until wrk is empty, any not-blacklisted neighbours within
        # cutoff distance from atoms that are blacklisted will be blacklisted
        # and pushed to wrk
        while len(wrk) > 0:
            atom_id = wrk.pop()
            
            # find neighbours approximately within cutoff
            ns = root.find_objects(atoms[rev_map[atom_id]].position,
                                     self._cutoff)
            for n in ns:
                if black_list[n]: # already blacklisted, ignore
                    continue
                if (self.length2(atoms[rev_map[n]].position,
                    atoms[rev_map[atom_id]].position) > cutoff2):
                    # outside cutoff distance, ignore
                    continue
                # not blacklisted, within cutoff, add to wrk
                # and blacklist
                wrk.append(n)
                black_list[n] = True

        k = 0
        subset = Atoms()
        # construct cleaned up subset of atoms within the CalcRegion
        for atom_id in root.get_objects():
            if black_list[atom_id]: # exclude blacklisted atom
                if self._debug > 2:
                    print("excluding atom: %i" % atom_id)
                continue
            subset.extend(atoms[rev_map[atom_id]:rev_map[atom_id]+1])
            weights.append(
                    self.get_weight(atoms[rev_map[atom_id]].position))
            subset_map[k] = rev_map[atom_id]
            k += 1

        # create weights array that matches subset
        wa = np.zeros((len(weights), 3))
        for i in range(len(weights)):
            wa[i] += weights[i]
        
        # set periodic boundary condition
        if self._pbc:
            subset.set_pbc(self._pbc)

        # dump CalcRegion contents to a trajectory file if debugging
        if self._debug_file != None:
            ase.io.pdb.write_pdb(self._debug_file, subset)

        return subset, subset_map, wa


class CalcBox(CalcRegion):
    """
    Atom selector for rectangular volumes with optional transition buffer
    region. The transition region adjusts weights starting from 0.0 on the
    outer surface and reaching 1.0 on the inner surface. The volume 
    contained within the inner surface gives weight 1.0.

    The volume can be made periodic by setting the pbc parameter for the
    constructor. It's simply passed on to the generated Atoms object, so it
    should adhere to the ASE documentation.
    """
    def __init__(self, name, pos=(0,0,0), dim=(1,1,1),
                 inner_dim=None, cutoff=2.0, pbc=None, debug=0):
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

        super(CalcBox, self).__init__(name, pos,
                            cutoff=cutoff, pbc=pbc, debug=debug)

        self._dim = dim

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


    def get_bounding_box(self):
        return (self._pos, self._dim)

    def get_weight(self, pos):
        x, y, z = pos

        # if outside the box, return 0.0
        if (x > self._xmax or x < self._xmin
                or y > self._ymax or y < self._ymin
                or z > self._zmax or z < self._zmin):
            return 0.0

        # if no inner box, straight to 1.0
        if not self._inner_dim:
            return 1.0
        
        # calculate transition region weight
        r = []

        if x >= self._xmin and x < self._inner_xmin:
            w = self._scale(x, self._xmin, self._inner_xmin)
            r.append(w)
        if x <= self._xmax and x > self._inner_xmax:
            w = self._scale(x, self._inner_xmax, self._xmax, reverse=1)
            r.append(w)
        if y >= self._ymin and y < self._inner_ymin:
            w = self._scale(y, self._ymin, self._inner_ymin)
            r.append(w)
        if y <= self._ymax and y > self._inner_ymax:
            w = self._scale(y, self._inner_ymax, self._ymax, reverse=1)
            r.append(w)
        if z >= self._zmin and z < self._inner_zmin:
            w = self._scale(z, self._zmin, self._inner_zmin)
            r.append(w)
        if z <= self._zmax and z > self._inner_zmax:
            w = self._scale(z, self._inner_zmax, self._zmax, reverse=1)
            r.append(w)

        if len(r) == 0:
            return 1.0

        weight = 1.0

        for w in r:
            weight = min(weight, w)

        return weight

    @staticmethod
    def _scale(value, left, right, reverse=0):
        if value < left:
            return 0.0
        if value > right:
            return 1.0
        if reverse == 0:
            return (value - left)/(right - left)
        else:
            return (right - value)/(right - left)


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


