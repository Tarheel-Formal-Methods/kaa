import sympy as sp
import numpy as np
from kaa.bundle import Bundle, BundleTransformer
from kaa.model import Model

def test_bund_trans_1():

    x,y = sp.Symbol('x'), sp.Symbol('y')
    dx = x + 1
    dy = y + 1
    dyns  =[dx, dy]
    vars = [x, y]

    L = np.empty([2,2])
    T = np.empty(2)

    L[0] = [1, 0]
    L[1] = [0, 1]

    T[0] = 0
    T[1] = 1

    T = [T]

    offu = np.empty(2)
    offl = np.empty(2)

    offu[0] = 1
    offu[1] = 1

    offl[0] = 1
    offl[1] = 1

    init_bund = Bundle(T, L, offu, offl, vars)
    trans = BundleTransformer(dyns)

    trans_bund = trans.transform(init_bund)

    print(trans_bund)
    assert False
