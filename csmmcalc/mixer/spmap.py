# Copyright (C) 2012 Aalto University
# Author: Lauri Leukkunen <lauri.leukkunen@aalto.fi>


class SpatialMap(object):
    def __init__(self, pos=(0,0,0), dim=(2,2,2), res=1.0):
        super(SpatialMap, self).__init__()
        self._pos = pos
        self._dim = dim
        self._submaps = None
        self._objects = []
        if (res < dim[0] or
            res < dim[1] or
            res < dim[2]):
            self._submaps = [[[0,0],[0,0]],[[0,0],[0,0]]]
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        new_x = pos[0] - (dim[0]/4.0) + i*dim[0]/2.0
                        new_y = pos[1] - (dim[1]/4.0) + j*dim[1]/2.0
                        new_z = pos[2] - (dim[2]/4.0) + k*dim[2]/2.0
                        self._submaps[i][j][k] = SpatialMap(
                            pos=(new_x, new_y, new_z),
                            dim=(dim[0]/2.0, dim[1]/2.0, dim[2]/2.0),
                            res=res)

    def add_object(self, obj, pos):
        if (pos[0] < self._pos[0] - self._dim[0]/2.0 or
            pos[0] > self._pos[0] + self._dim[0]/2.0 or
            pos[1] < self._pos[1] - self._dim[1]/2.0 or
            pos[1] > self._pos[1] + self._dim[1]/2.0 or
            pos[2] < self._pos[2] - self._dim[2]/2.0 or
            pos[2] > self._pos[2] + self._dim[2]/2.0):
            return
        if not self._submaps:
            self._objects.append(obj)
        else:
            self.get_submap(pos).add_object(obj, pos)
    
    def get_submap(self, pos):
        if not self._submaps:
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
        return self._submaps[i][j][k]

    def get_objects(self):
        if not self._submaps:
            return self._objects
        res = []
        for s1 in self._submaps:
            for s2 in s1:
                for s3 in s2:
                    res += s3.get_objects()
        return res

    def _get_flags(self, pos, coord, reach):
        flags = [0, 0]
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
        if not self._submaps:
            return self._objects
        
        x_flags = self._get_flags(pos, 0, reach)
        y_flags = self._get_flags(pos, 1, reach)
        z_flags = self._get_flags(pos, 2, reach)

        res = [] 
   
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    if x_flags[i] and y_flags[j] and z_flags[k]:
                        res += (
                            self._submaps[i][j][k].find_objects(pos, reach))
        return res


if __name__ == "__main__":
    edge = 100.0
    count = 100
    sm = SpatialMap((0,0,0), (edge, edge, edge), res=10.0)
    import numpy as np
    x = np.linspace(-edge/2.0, edge/2.0, count)
    y = np.linspace(-edge/2.0, edge/2.0, count)
    z = np.linspace(-edge/2.0, edge/2.0, count)

    print("spatial map created")
    for xi in x:
        for yi in y:
            for zi in z:
                pos = (xi, yi, zi)
                sm.add_object("object: (%f, %f, %f)" % pos, pos)

    print("objects added")
    print("total objects: %i" % len(sm.get_objects()))
    subset = sm.find_objects((.5,.5,.5), 0.33)
    print("found objects: %i" % len(subset))
