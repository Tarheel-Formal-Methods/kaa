import sympy as sp
import numpy as np

from kaa.bundle import Bundle
from kaa.model import Model

class Phosphorelay(Model):

    def __init__(self):

        dim_sys = 7

        x1, x2, x3, x4, x5, x6, x7 = sp.Symbol("x1"), sp.Symbol("x2"), sp.Symbol("x3"), sp.Symbol("x4"), sp.Symbol("x5"), sp.Symbol("x6"), sp.Symbol("x7");
        vars = [x1, x2, x3, x4, x5, x6, x7]

        delta = 0.01

        dx1 = x1 + ( -0.4*x1 + 5*x3*x4 )*delta
        dx2 = x2 + ( 0.4*x1 - x2 )*delta
        dx3 = x3 + ( x2-5*x3*x4 )*delta
        dx4 = x4 + ( 5*x5*x6 - 5*x3*x4 )*delta
        dx5 = x5 + ( -5*x5*x6 + 5*x3*x4 )*delta
        dx6 = x6 + ( 0.5*x7 - 5*x5*x6 )*delta
        dx7 = x7 + ( -0.5*x7 + 5*x5*x6 )*delta

        dyns = [dx1,dx2,dx3,dx4,dx5,dx6,dx7];

        num_dirs = 10
        num_temps = 2

        L = np.zeros([num_dirs,dim_sys])
        for i in range(dim_sys):
            L[i][i] = 1

        L[7][2] = 1; L[7][3] = 1;
        L[8][4] = 1; L[8][5] = 1;
        L[9][2] = 1; L[9][3] = 1; L[9][4] = 1; L[9][5] = 1;

        T = np.zeros([num_temps,dim_sys])
        T[0][0] = 0; T[0][1] = 1; T[0][2] = 2; T[0][3] = 3; T[0][4] = 4; T[0][5] = 5; T[0][6] = 6;
        T[1][0] = 0; T[1][1] = 1; T[1][2] = 2; T[1][3] = 7; T[1][4] = 4; T[1][5] = 5; T[1][6] = 6;
        #T[2][0] = 0; T[2][1] = 1; T[2][2] = 2; T[2][3] = 7; T[2][4] = 8; T[2][5] = 5; T[2][6] = 6;
        #T[3][0] = 0; T[3][1] = 1; T[3][2] = 2; T[3][3] = 7; T[3][4] = 9; T[3][5] = 5; T[3][6] = 6;

        offu = np.zeros(num_dirs)
        offl = np.zeros(num_dirs)

        offu[0] = 1.01; offl[0] = -1.00;
        offu[1] = 1.01; offl[1] = -1.00;
        offu[2] = 1.01; offl[2] = -1.00;
        offu[3] = 1.01; offl[3] = -1.00;
        offu[4] = 1.01; offl[4] = -1.00;
        offu[5] = 1.01; offl[5] = -1.00;
        offu[6] = 1.01; offl[6] = -1.00;

        offu[7] = 100; offl[7] = 100;
        offu[8] = 100; offl[8] = 100;
        offu[9] = 100; offl[9] = 100;


        B = Bundle(T,L,offu,offl,vars)

        super().__init__(B, dyns, vars)
