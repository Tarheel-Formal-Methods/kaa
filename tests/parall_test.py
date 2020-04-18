import numpy as np
import sympy as sp
from sapo.parallelotope import Parallelotope

def test_min_max_parall_1():

    x, y = sp.Symbol('x'), sp.Symbol('y')

    A = np.empty([4,2])
    b = np.empty(4)

    A[0] = [1,0]
    A[1] = [0,1]
    A[2] = [-1,0]
    A[3] = [0,-1]

    b[0] = 1
    b[1] = 1
    b[2] = 1
    b[3] = 1

    p = Parallelotope(A, b, [x,y])
    min_point, max_point = p.getMinPoint(), p.getMaxPoint()

    assert np.equal(max_point, [1,1]).all()
    assert np.equal(min_point, [-1,-1]).all()
