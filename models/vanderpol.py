import sympy as sp
import numpy as np

from kaa.bundle import Bundle
from kaa.model import Model


class VanDerPol(Model):

    def __init__(self):

        x1, y1, x2, y2 = sp.Symbol('x1'), sp.Symbol('y1'),  sp.Symbol('x2'), sp.Symbol('y2')
        delta = 0.1
        dx1 =  x1 + y1*delta
        dy1 = y1 + (1*(1-x1**2)*y1 - 2*x1 + x2)*delta
        dx2 = x2 + y2*delta
        dy2 = y2 + (2*(1 - x2**2)*y2 - 2*x2 + x1)*delta

        vars = [x1,y1, x2, y2]
        dyns = [dx1, dy1, dx1, dy2]

        dim_sys = 4
        num_dirs = 4
        num_temps = 3

        L = np.zeros([num_dirs, dim_sys])

        L[0][0] = 1;
        L[1][1] = 1;
        L[2][2] = 1;
        L[3][3] = 1;

        T = np.zeros([num_temps, dim_sys])
        T[0][0] = 0; T[0][1] = 1;
        T[1][0] = 1; T[1][1] = 2;
        T[2][0] = 2; T[2][1] = 3;
        #T[5][0] = 2; T[5][1] = 3;

        offu = np.zeros(num_dirs);
        offl = np.zeros(num_dirs);

        offu[0] = 1.85; offl[0] = -1.55;
        offu[1] = 2.45; offl[1] = -2.35;
        offu[2] = 1.85; offl[2] = -1.55;
        offu[3] = 2.45; offl[3] = -2.35;

        b = Bundle(T, L,  offu, offl, vars)
        super().__init__(b, dyns, vars)
