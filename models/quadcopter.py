import sympy as sp
import numpy as np

from kaa.bundle import Bundle
from kaa.model import Model


class Quadcopter(Model):

    def __init__(self):

           dim_sys = 17
           pn, pe, h, u, v, w, q0v, q1v = sp.Symbol("pn"), sp.Symbol("pe"), sp.Symbol("h"), sp.Symbol("u"), sp.Symbol("v"), sp.Symbol("w"), sp.Symbol("q0v"), sp.Symbol("q1v")
           q2v, q3v, p, q, r, hI, uI, vI, psiI = sp.Symbol("q2v"), sp.Symbol("q3v"), sp.Symbol("p"),sp.Symbol("q"),sp.Symbol("r"), sp.Symbol("hI"), sp.Symbol("uI"), sp.Symbol("vI"), sp.Symbol("psiI")

           vars = [pn, pe, h, u, v, w, q0v, q1v, q2v, q3v, p, q, r, hI, uI, vI, psiI]

           M = 0.0015;
           mr = 0.001;
           R = 0.020;
           l = 0.045;
           g = 9.81;
           m = M + 4*mr;
           Jx = (2*M*R**2)/5 + (2*l**2)*mr;
           Jy = (2*M*R**2)/5 + (2*l**2)*mr;
           Jz = (2*M*R**2)/5 + (4*l**2)*mr;

           ur = 0;
           vr = 0;
           psir = 0;
           hr = 1;

           phi = 2*q1v;
           theta = 2*q2v;
           psi = 2*q3v;

           delta = 0.01;

           F = 0.0361*hI + 0.0694*h + 0.0603*w;
           tauphi = -0.0003*vI - 0.0005*v - 0.0018*phi - 0.0004*p;
           tautheta = 0.0003*uI + 0.0005*u - 0.0018*theta - 0.0004*q;
           taupsi = -0.0003*psiI - 0.0006*psi - 0.0003*r;

           dpn = pn + (u*(2*q0v**2 + 2*q1v**2 - 1) - v*(2*q0v*q3v - 2*q1v*q2v ) + w*(2*q0v*q2v + 2*q1v*q3v ))*delta;
           dpe = pe + (v*(2*q0v**2 + 2*(q2v**2) - 1) + u*(2*q0v*q3v + 2*q1v*q2v ) - w*(2*q0v*q1v - 2*q2v*q3v ))*delta;
           dh = h + (w*(2*(q0v**2) + 2*(q3v**2) - 1) - u*(2*q0v*q2v - 2*q1v*q3v ) + v*(2*q0v*q1v + 2*q2v*q3v ))*delta;

           du = u + (r*v - q*w - g*(2*q0v*q2v - 2*q1v*q3v ))*delta;
           dv = v + (p*w - r*u + g*(2*q0v*q1v + 2*q2v*q3v ))*delta;
           dw = w + (q*u - p*v - F/m + g*(2*(q0v**2) + 2*(q3v**2) - 1 ))*delta;

           dq0v = q0v +(-(q1v/2)*p - (q2v/2)*q - (q3v/2)*r)*delta;
           dq1v = q1v + ((q0v/2)*p - (q3v/2)*q + (q2v/2)*r)*delta;
           dq2v = q2v + ((q3v/2)*p + (q0v/2)*q - (q1v/2)*r)*delta;
           dq3v = q3v + ((q1v/2)*q - (q2v/2)*p + (q0v/2)*r)*delta;

           dp = p + ((1/Jx)*tauphi + ((Jy - Jz)/Jx)*q*r)*delta;
           dq = q + ((1/Jy)*tautheta - ((Jx - Jz)/Jy)*p*r)*delta;
           dr = r + ((1/Jz)*taupsi + ((Jx - Jy)/Jz)*p*q)*delta;

           dhI = hI + (h - hr)*delta;
           duI = uI +(u - ur)*delta;
           dvI = vI + (v - vr)*delta;
           dpsiI = psiI + (psi - psir)*delta;

           dyns = [dpn,dpe,dh,du,dv,dw,dq0v,dq1v,dq2v,dq3v,dp,dq,dr,dhI,duI,dvI,dpsiI]

           num_dirs = 18
           num_temp = 2

           L = np.zeros([num_dirs,dim_sys])

           for i in range(dim_sys):
               L[i][i] = 1

           L[17][2] = 0.5; L[17][5] = 0.5; L[17][6] = 0.5; L[17][15] = 0.25;

           T = np.zeros([num_temp, dim_sys]);
           for i in range(dim_sys):
               T[0][i] = i
               T[1][i] = i

           T[1][5] = 17

           offu = np.zeros(num_dirs);
           offl = np.zeros(num_dirs);

           offu[2] = 0.21; offl[2] = -0.20;
           offu[6] = 1; offl[6] = -1;
           offu[17] = 100; offl[17] = 100;

           b = Bundle(T, L, offu, offl, vars)
           super().__init__(b, dyns, vars)
