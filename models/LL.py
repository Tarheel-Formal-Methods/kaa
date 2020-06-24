import sympy as sp
import numpy as np

from kaa.bundle import Bundle
from kaa.model import Model


class LL(Model):

    def __init__(self):

        x1, x2, x3, x4, x5, x6, x7 = sp.Symbol('x1'), sp.Symbol('x2'), sp.Symbol('x3'), sp.Symbol('x4'),sp.Symbol('x5'), sp.Symbol('x6'), sp.Symbol('x7')
        delta = 0.05

        dx1 =  x1 + (1.4*x3 - 0.9*x1)*delta
        dx2 = x2 + (2.5*x5 - 1.5*x2)*delta
        dx3 = x3 + (0.6*x7 - 0.8*x2*x3)*delta
        dx4 = x4 + (2 - 1.3*x3*x4)*delta
        dx5 = x5 + (0.7*x1 - x4*x5)*delta
        dx6 = x6 + (0.3*x1 - 3.1*x6)*delta
        dx7 = x7 + (1.8*x6 - 1.5*x2*x7)*delta

        vars = [x1, x2, x3, x4, x5, x6, x7]
        dyns = [dx1, dx2, dx3, dx4, dx5, dx6, dx7]

        dim_sys = 7
        num_dirs = 9
        num_temps = 2

        L = np.zeros([num_dirs, dim_sys])

        for i in range(dim_sys):
            L[i][i] = 1

        L[7][2] = 1; L[7][3] = 1;
        L[8][3] = 1; L[8][4] = 1;


        T = np.zeros([num_temps, dim_sys])
        T[0][0] = 0; T[0][1] = 1; T[0][2] = 2; T[0][3] = 3; T[0][4] = 4; T[0][5] = 5; T[0][6] = 6;
        #T[0][0] = 1; T[0][1] = 2; T[0][2] = 3; T[0][3] = 4; T[0][4] = 5; T[0][5] = 6; T[0][6] = 7;
        T[1][0] = 2; T[1][1] = 3; T[1][2] = 4; T[1][3] = 5; T[1][4] = 6; T[1][5] = 7; T[1][6] = 8;

        offu = np.zeros(num_dirs)
        offl = np.zeros(num_dirs)

        '''
        offu[0] = 1.3; offl[0] = -1.1;
        offu[1] = 1.15; offl[1] = -0.95;
        offu[2] = 1.6; offl[2] = -1.4;
        offu[3] = 2.5; offl[3] = -2.3;
        offu[4] = 1.1; offl[4] = -0.9;
        offu[5] = 0.2; offl[5] = 0.0;
        offu[6] = 0.55; offl[6] = -0.35;
        '''
        offu[0] = 1.21; offl[0] = -1.19;
        offu[1] = 1.06; offl[1] = -1.04;
        offu[2] = 1.51; offl[2] = -1.49;
        offu[3] = 2.41; offl[3] = -2.39;
        offu[4] = 1.01; offl[4] = -0.99;
        offu[5] = 0.11; offl[5] = -0.09;
        offu[6] = 0.46; offl[6] = -0.44;
        offu[7] = 2.87; offl[7] = -2.81;
        offu[8] = 2.66; offl[8] = -2.58;

        b = Bundle(T, L,  offu, offl, vars)
        super().__init__(b, dyns, vars)
