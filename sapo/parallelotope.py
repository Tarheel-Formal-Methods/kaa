import numpy as np
from scipy.optimize import linprog


class Parallelotope:

    def __init__(self, A, b, vars):
        self.A = A
        self.b = b
        self.vars = vars

    def getMinPoint(self):

        c = [1 for _ in enumerate(self.vars)]
        var_min = linprog(c, A_ub=self.A, b_ub=self.b)

        return var_min.x

    def getMaxPoint(self):

        c = [-1 for _ in enumerate(self.vars)]
        var_max = linprog(c, A_ub=self.A, b_ub=self.b)

        return var_max.x
