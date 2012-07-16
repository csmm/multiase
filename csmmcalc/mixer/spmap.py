# Copyright (C) 2012 Aalto University
# Author: Lauri Leukkunen <lauri.leukkunen@aalto.fi>


class SpatialMap(object):
    def __init__(self, pos=(0,0,0), dim=(2,2,2), res=1.0):
        super(SpatialMap, self).__init__()
        self._pos = pos
        self._dim = dim

        

