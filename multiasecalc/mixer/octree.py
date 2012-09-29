# Copyright (C) 2012 Aalto University
# Author: Lauri Leukkunen <lauri.leukkunen@aalto.fi>


class OctreeNode(object):
    """
    OctreeNode divides 3D cartesian space into 2x2x2 grid of 8 blocks which are
    bound to child nodes (subnodes in this source file). Child nodes are only
    created as needed and each object stored in the datastructure will reside in
    a leaf node. The degree (or depth) of the Octree is determined by the
    resolution desired.  More accurate octrees require a deeper chain of inner
    nodes to reach a granularity of space that is equal or better than the given
    resolution.

    Each node's position in space is given by the central point of the volume it
    occupies. The dimensions are given as the full length of each side of the 3D
    box.

    Objects are stored in python set()s, so each leaf node can contain the same
    object only once. This is done for performance reasons.
    """
    def __init__(self, pos=(0,0,0), dim=(2,2,2), res=1.0):
        super(OctreeNode, self).__init__()
        self._pos = pos
        self._dim = dim
        self._res = res
        self._subnodes = None
        self._objects = set()
        if (res < dim[0] or
            res < dim[1] or
            res < dim[2]):
            self._subnodes = [[[0,0],[0,0]],[[0,0],[0,0]]]


    def _create_subnode(self, addr):
        i, j, k = addr
        new_x = self._pos[0] - (self._dim[0]/4.0) + i*self._dim[0]/2.0
        new_y = self._pos[1] - (self._dim[1]/4.0) + j*self._dim[1]/2.0
        new_z = self._pos[2] - (self._dim[2]/4.0) + k*self._dim[2]/2.0
        return OctreeNode(
                pos=(new_x, new_y, new_z),
                dim=(self._dim[0]/2.0, self._dim[1]/2.0, self._dim[2]/2.0),
                res=self._res)


    def add_object(self, obj, pos):
        # check pos is within the node
        if (pos[0] < self._pos[0] - self._dim[0]/2.0 or
            pos[0] > self._pos[0] + self._dim[0]/2.0 or
            pos[1] < self._pos[1] - self._dim[1]/2.0 or
            pos[1] > self._pos[1] + self._dim[1]/2.0 or
            pos[2] < self._pos[2] - self._dim[2]/2.0 or
            pos[2] > self._pos[2] + self._dim[2]/2.0):
            return
        if self._subnodes == None:
            self._objects.add(obj)
            return

        i,j,k = self.get_subnode_addr(pos)
        if self._subnodes[i][j][k] == 0:
            self._subnodes[i][j][k] = self._create_subnode((i,j,k))
        
        self._subnodes[i][j][k].add_object(obj, pos)
    
    def get_subnode_addr(self, pos):
        if not self._subnodes:
            return None

        i = 0
        j = 0
        k = 0
        if pos[0] > self._pos[0]:
            i = 1
        if pos[1] > self._pos[1]:
            j = 1
        if pos[2] > self._pos[2]:
            k = 1
        return (i,j,k)

    def get_objects(self):
        if not self._subnodes:
            return self._objects
        res = set()
        for s1 in self._subnodes:
            for s2 in s1:
                for s3 in s2:
                    if s3 == 0:
                        continue
                    res = res.union(s3.get_objects())
        return res

    def _get_flags(self, pos, coord, reach):
        flags = [0, 0]
        
        # first check if search is within our volume
        if pos[coord] + reach < self._pos[coord] - self._dim[coord]/2.0:
            return flags
        if pos[coord] - reach > self._pos[coord] + self._dim[coord]/2.0:
            return flags

        if pos[coord] < self._pos[coord]:
            if pos[coord] + reach > self._pos[coord]:
                # get both splits
                flags[0] = 1
                flags[1] = 1
            else:
                # just lower
                flags[0] = 1
        else:
            if pos[coord] - reach < self._pos[coord]:
                # get both splits
                flags[0] = 1
                flags[1] = 1
            else:
                # just higher
                flags[1] = 1
        return flags


    def find_objects(self, pos, reach=1.0):
        if not self._subnodes:
            return self._objects
        
        x_flags = self._get_flags(pos, 0, reach)
        y_flags = self._get_flags(pos, 1, reach)
        z_flags = self._get_flags(pos, 2, reach)

        res = set()
   
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    if x_flags[i] and y_flags[j] and z_flags[k]:
                        if self._subnodes[i][j][k] == 0:
                            continue
                        res = res.union(
                            self._subnodes[i][j][k].find_objects(pos, reach))
        return res


