import numpy as np
import scipy as sp


class Parallelotope:

    def __init__(self, A, b, vars):
        self.A = A
        self.b = b
        self.vars = vars

    def getMinPoint(self):

        min_coord = []
        for var_ind, _ in enumerate(self.vars):
            c = [0 for _ in enumerate(self.vars)]
            c[var_ind] = 1
            var_min = sp.optimize.linprog(c,self.A,self.b)
            min_coord.append(var_min.fun)

        return min_coord

    def getMaxPoint(self):

        max_coord = []
        for var_ind, _ in enumerate(self.vars):
            c = [0 for _ in enumerate(self.vars)]
            c[var_ind] = -1
            var_max = sp.optimize.linprog(c,self.A,self.b)
            max_coord.append(-1 * var_max.fun)

        return max_coord
