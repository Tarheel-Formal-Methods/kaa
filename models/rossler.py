import sympy as sp
import numpy as np

from kaa.model import Model
from kaa.bundle import Bundle

class Rossler(Model):

    def __init__(self):

        x, y, z = sp.Symbol('x'), sp.Symbol('y'), sp.Symbol('z')
        vars = [x, y, z]

        dim_sys = len(vars)

        dx = x + (-y-z)*0.025
        dy = y + (x + 0.1*y)*0.025
        dz = z + (0.1 + z*(x-14))*0.02

        dyns = [dx, dy ,dz]

        num_direct = 5
        num_temps = 3

        L = np.zeros([num_direct, dim_sys])
        T = np.zeros([num_temps, dim_sys])

        L[0][0] = 1
        L[1][1] = 1
        L[2][2] = 1

        L[3][0] = 1
        L[3][1] = 0.5

        L[4][0] = 0.5
        L[4][2] = 0.5

        T[0][0] = 0; T[0][1] = 1; T[0][2] = 2;
        T[1][0] = 1; T[1][1] = 2; T[1][2] = 3;
        T[2][0] = 2; T[2][1] = 3; T[2][2] = 4;

        offu = np.zeros(num_direct)
        offl = np.zeros(num_direct)

        offu[0] = 0.1; offl[0] = -0.09;
        offu[1] = 5; offl[1] = -4.99;
        offu[2] = 0.1; offl[2] = -0.09;
        offu[3] = 10; offl[3] = 0;
        offu[4] = 10; offl[4] = 0;

        b = Bundle(T, L, offu, offl, vars)

        super().__init__(b, dyns, vars)
