import numpy as np
from kaa.lputil import minLinProg, maxLinProg
"""
Object encapsulating routines calculating properties of parallelotopes
"""
class Parallelotope:

    def __init__(self, A, b, vars):
        self.A = A
        self.b = b
        self.vars = vars

    '''
     Returns the minimum base point of the parallelotope.
    '''
    def getMinPoint(self):

        var_min = []
        for var_ind,_ in enumerate(self.vars):
            c = [0 for _ in enumerate(self.vars)]
            c[var_ind] = 1
            min_point = minLinProg(c, self.A, self.b).fun
            var_min.append(min_point)
        return var_min

    '''
    Returns the maximum base point of the parallelotope.
    '''
    def getMaxPoint(self):

        var_max = []
        for var_ind,_ in enumerate(self.vars):
            c = [0 for _ in enumerate(self.vars)]
            c[var_ind] = 1
            max_point = maxLinProg(c, self.A, self.b).fun
            var_max.append(max_point)
        return var_max
