import numpy as np
from sapo.lputil import minLinProg, maxLinProg


class Parallelotope:

    def __init__(self, A, b, vars):
        self.A = A
        self.b = b
        self.vars = vars
        self.c = [1 for _ in enumerate(self.vars)]

    def getMinPoint(self):

        var_min = minLinProg(self.c, self.A, self.b)
        return var_min.x

    def getMaxPoint(self):

        var_max = maxLinProg(self.c, self.A, self.b)
        return var_max.x
