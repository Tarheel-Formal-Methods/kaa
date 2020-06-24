import sympy as sp
import numpy as np

from kaa.bundle import Bundle
from kaa.model import Model

class Ebola(Model):

    def __init__(self):
        dim_sys = 5

        s, e, q, i, r, kappa1, gamma1 =  sp.Symbol("s"), sp.Symbol("e"), sp.Symbol("q"), sp.Symbol("i"), sp.Symbol("r"), sp.Symbol("kappa1"), sp.Symbol("gamma1");
        vars = [s, e,q, i, r]
        params = [kappa1, gamma1]

        beta = 0.35;
        kappa2 = 0.3;
        gamma2 = 0.6;
        sigma = 0.28;
        Delta = 0.5;

        ds = s - ((s*beta*i) + gamma1*q)*Delta;
        de = e + ((s*beta*i) - (kappa1+kappa2)*e)*Delta;
        dq = q + (kappa1*e - (gamma1+gamma2)*q)*Delta;
        di = i + (gamma2*q + kappa2*e - sigma*i)*Delta;
        dr = r + (sigma*i)*Delta;

        dyns = [ds,de,dq,di,dr]

        num_dirs = 5

        L = np.zeros([num_dirs, dim_sys])
        L[0][0] = 1;
        L[1][1] = 1;
        L[2][2] = 1;
        L[3][3] = 1;
        L[4][4] = 1;

        offu = np.zeros(num_dirs)
        offl = np.zeros(num_dirs)

        offu[0] = 0.8; offl[0] = -0.79;
        offu[1] = 0; offl[1] = 0;
        offu[2] = 0.0000; offl[2] = 0.00000;
        offu[3] = 0.2; offl[3] = -0.19;
        offu[4] = 0; offl[4] = 0;

        T = np.zeros([1,dim_sys])

        T[0][0] = 0; T[0][1] = 1; T[0][2] = 2; T[0][3] = 3; T[0][4] = 4;

        B = Bundle(T, L, offu, offl, vars)

        super().__init__(B, dyns, vars)
