# Copyright (C) 2012 Aalto University
# Author: Lauri Leukkunen <lauri.leukkunen@aalto.fi>

from csmmcalc.mixer.octree import OctreeNode
from csmmcalc.mixer.mixer import Mixer
from csmmcalc.mixer.selector import CalcBox
from ase.atoms import Atoms
from ase.structure import molecule


import unittest
import numpy as np
import math


class IntContainer(object):
    """
    Utility class to help with OctreeNode testing
    """
    pos = (0.0, 0.0, 0.0)
    value = 0

    def distance2(self, a):
        dx = self.pos[0] - a.pos[0]
        dy = self.pos[1] - a.pos[1]
        dz = self.pos[2] - a.pos[2]
        return dx**2 + dy**2 + dz**2

    def distance(self, a):
        return math.sqrt(self.distance2(b))


origo = IntContainer()

class OctreeTests(unittest.TestCase):
    def setUp(self):
        self._edge = 100.0
        self._count = 5
        self._res = 5.0
        self._root = OctreeNode((0,0,0),
                (self._edge, self._edge, self._edge),
                res=5.0)
        
        x = np.linspace(-self._edge/2.0, self._edge/2.0, self._count)
        
        c = 0
        for i in x:
            for j in x:
                for k in x:
                    tmp = IntContainer()
                    tmp.pos = (i,j,k)
                    tmp.value = c
                    self._root.add_object(tmp, tmp.pos)
                    #print("%i: (%f,%f,%f)" % (c, i, j, k))
                    c += 1

    def test_octree_find_objects_1(self):
        correct_values = [62]
        subset = self._root.find_objects((0., 0., 0.), 10.)
        values = [x.value for x in subset]
        self.assertEqual(set(values), set(correct_values))

    def test_octree_find_objects_2(self):
        correct_values = [31, 32, 33, 36, 37, 38, 41, 42,
                43, 56, 57, 58, 61, 62, 63, 66, 67, 68, 81,
                82, 83, 86, 87, 88, 91, 92, 93]
        subset = self._root.find_objects((0., 0., 0.), 26.)
        values = [x.value for x in subset]
        self.assertEqual(set(values), set(correct_values))


class SelectorTests(unittest.TestCase):

    def _set_atoms_pos(self, atoms, pos):
        p = atoms.get_positions()
        p = p + np.array(pos)
        atoms.set_positions(p)

    def _set_molecule_pos(self, atoms, molecule, pos):
        p = atoms.get_positions()
        dx = pos[0] - p[molecule[0]][0]
        dy = pos[1] - p[molecule[0]][1]
        dz = pos[2] - p[molecule[0]][2]
        for i in molecule:
            p[i][0] += dx
            p[i][1] += dy
            p[i][2] += dz
        atoms.set_positions(p)

    def setUp(self):
        self._atoms = Atoms()
        self._mcount = 5
        self._edge = 100.0

        x = np.linspace(-self._edge/2.0, self._edge/2.0, self._mcount)
        c = 0
        for i in x:
            for j in x:
                for k in x:
                    ch4 = molecule("CH4")
                    self._set_atoms_pos(ch4, (i,j,k))
                    self._atoms.extend(ch4)
                    c += 1
        
        Mixer.set_atom_ids(self._atoms)


    def test_calcbox_all(self):
        cb = CalcBox(pos=(0., 0., 0.),
                     dim=(self._edge + 5, self._edge + 5, self._edge + 5))
        sa, srmap, sw = cb.select_atoms(self._atoms)
        ra = self._atoms
        self.assertEqual(ra, sa)

    def test_calcbox_subset(self):
        cb = CalcBox(pos=(0., 0., 0.),
                     cutoff=0.0,
                     dim=(50.0, 50.0, 50.0))
        sa, srmap, sw = cb.select_atoms(self._atoms)
        correct_ids = set([307, 215, 405, 280, 281, 410, 155, 156, 285, 286,
            415, 160, 289, 290, 419, 165, 294, 455, 430, 432, 305, 306, 435,
            180, 414, 310, 439, 312, 185, 186, 315, 188, 190, 181, 193, 343,
            161, 457, 330, 311, 332, 205, 462, 335, 337, 210, 340, 213, 313,
            218, 314, 465, 444, 338, 437, 440, 318, 319, 460])
        self.assertEqual(set(Mixer.get_atom_ids(sa)), correct_ids)

    def test_calcbox_moving_ch4(self):
        cb = CalcBox(pos=(0., 0., 0.),
                     cutoff=2.0,
                     dim=(55.0, 55.0, 55.0),
                     debug=2)
        ch4 = [0, 1, 2, 3] # these are atoms indexes and ids
        z = np.linspace(-self._edge/2.0 - 25, self._edge/2.0 + 25, 100)
        x = np.zeros_like(z)
        y = np.zeros_like(z)
        pos = np.column_stack((x, y, z))
        for i in range(len(pos)):
            print(pos[i])
            self._set_molecule_pos(self._atoms,
                                   ch4,
                                   pos[i])
            sa, srmap, sw = cb.select_atoms(self._atoms)
            s_ids = set(Mixer.get_atom_ids(sa))
            all_ids = set(Mixer.get_atom_ids(self._atoms))
            inside_ids = s_ids.intersection(set(ch4))
            if not (len(inside_ids) == 0 or len(inside_ids) == 4):
                print("CH4 pos: %s" % pos[i])
                print("CH4 in selection: %s" % inside_ids)

    def test_weights(self):
        cb = CalcBox(pos=(0., 0., 0.),
                     cutoff=2.0,
                     dim=(100.0, 100.0, 100.0),
                     inner_dim=(50.0, 50.0, 50.0))
        # test all three axis
        test_positions = np.column_stack((np.linspace(-55.0, 55.0, 100),
                                          np.zeros(100),
                                          np.zeros(100)))
        test_positions = np.append(test_positions, np.column_stack(
                                (np.zeros(100),
                                 np.linspace(-55.0, 55.0, 100),
                                 np.zeros(100))), axis=0)
        test_positions = np.append(test_positions, np.column_stack(
                                (np.zeros(100),
                                 np.zeros(100),
                                 np.linspace(-55.0, 55.0, 100))),
                                 axis=0)
        # and all 8 corners too
        a = 37.5
        test_positions = np.append(test_positions, np.array((
                                (-a, -a, -a),
                                (-a, -a, a),
                                (-a, a, -a),
                                (-a, a, -a),
                                (a, -a, -a),
                                (a, -a, a),
                                (a, a, -a),
                                (a, a, a))),
                                axis=0)
        correct_weights = np.loadtxt("weight_test.csv", delimiter=",")
        for i in range(len(test_positions)):
            w = cb.get_weight(test_positions[i])
            self.assertEqual(correct_weights[i], w)


if __name__ == "__main__":
    unittest.main()

