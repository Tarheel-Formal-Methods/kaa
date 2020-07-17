import sympy as sp
import numpy as np

from kaa.bundle import Bundle
from kaa.model import Model

class Basic2(Model):

    def __init__(self):

        x,y = sp.Symbol('x'), sp.Symbol('y')

        dx = x - 1
        dy = y - 1

        dyns  = [dx, dy]
        vars = [x, y]

        L = np.empty([2,2])
        T = np.empty([1,2])

        L[0] = [1, 0]
        L[1] = [0, 1]

        T[0][0] = 0
        T[0][1] = 1

        offu = np.empty(2)
        offl = np.empty(2)

        offu[0] = 1
        offu[1] = 1

        offl[0] = 1
        offl[1] = 1
    
        b = Bundle(T, L, offu, offl, vars)

        super().__init__(b, dyns, vars)
