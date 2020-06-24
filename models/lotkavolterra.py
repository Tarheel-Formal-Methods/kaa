import sympy as sp
import numpy as np

from kaa.bundle import Bundle
from kaa.model import Model


class LotkaVolterra(Model):

    def __init__(self):
        dim_sys = 5
        x1, x2, x3, x4, x5 = sp.Symbol("x1"), sp.Symbol("x2"), sp.Symbol("x3"), sp.Symbol("x4"), sp.Symbol("x5")
        vars = [x1, x2, x3, x4, x5]

        alpha = 0.85
        beta = 0.5
        delta = 0.01


        dx1 = x1 + ( x1*(1 - (x1 + alpha*x2 + beta*x5)) )*delta
        dx2 = x2 + (x2*(1 - (x2 + alpha*x3 + beta*x1)) )*delta
        dx3 = x3 + (x3*(1 - (x3 + alpha*x4 + beta*x2)) )*delta
        dx4 = x4 + (x4*(1 - (x4 + alpha*x5 + beta*x3)) )*delta
        dx5 = x5 + (x5*(1 - (x5 + alpha*x1 + beta*x4)) )*delta

        dyns = [dx1, dx2, dx3, dx4, dx5]

        num_dirs = 7
        num_temps = 1

        L = np.zeros([num_dirs,dim_sys])
        T = np.zeros([num_temps,dim_sys])

        for i in range(dim_sys):
            L[i][i] = 1

        L[5][0] = 1; L[5][1] = 1; L[5][2] = 1;
        L[6][3] = -1; L[6][4] = 1; L[6][0] = -1;

        T[0][0] = 0; T[0][1] = 1; T[0][2] = 2; T[0][3] = 3; T[0][4] = 4;
        #T[1][0] = 5; T[1][1] = 6; T[1][2] = 1; T[1][3] = 2; T[1][4] = 3;
        #T[2][0] = 5; T[2][1] = 6; T[2][2] = 2; T[2][3] = 3; T[2][4] = 4;

        offu = np.zeros(num_dirs)
        offl = np.zeros(num_dirs)

        offu[0] = 1.0; offl[0] = -0.95;
        offu[1] = 1.0; offl[1] = -0.95;
        offu[2] = 1.0; offl[2] = -0.95;
        offu[3] = 1.0; offl[3] = -0.95;
        offu[4] = 1.0; offl[4] = -0.95;

        offu[5] = 10.0; offl[5] = 1;
        offu[6] = 10.0; offl[6] = 1;

        b = Bundle(T,L,offu,offl,vars)

        super().__init__(b, dyns, vars)
