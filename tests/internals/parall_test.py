import numpy as np
import sympy as sp
from kaa.parallelotope import Parallelotope

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


def test_min_max_parall_2():

    x, y = sp.Symbol('x'), sp.Symbol('y')

    A = np.empty([4,2])
    b = np.empty(4)

    A[0] = [1,1]
    A[1] = [-1,1]
    A[2] = [-1,-1]
    A[3] = [1,-1]

    b[0] = 2
    b[1] = 2
    b[2] = 2
    b[3] = 2

    p = Parallelotope(A, b, [x,y])
    min_point, max_point = p.getMinPoint(), p.getMaxPoint()

    assert np.equal(max_point, [1,1]).all()
    assert np.equal(min_point, [-1,-1]).all()

def test_min_max_parall_3():

    x, y = sp.Symbol('x'), sp.Symbol('y')

    A = np.empty([4,2])
    b = np.empty(4)

    A[0] = [1,0]
    A[1] = [0,1]
    A[2] = [-1,0]
    A[3] = [0,-1]

    b[0] = 0
    b[1] = 0
    b[2] = 2
    b[3] = 2

    p = Parallelotope(A, b, [x,y])
    min_point, max_point = p.getMinPoint(), p.getMaxPoint()

    assert np.equal(max_point, [0,0]).all()
    assert np.equal(min_point, [-2,-2]).all()

def test_min_max_parall_4():

    x, y = sp.Symbol('x'), sp.Symbol('y')

    A = np.empty([4,2])
    b = np.empty(4)

    A[0] = [1,0]
    A[1] = [0,1]
    A[2] = [-1,0]
    A[3] = [0,-1]

    b[0] = -2
    b[1] = -2
    b[2] = 4
    b[3] = 4

    p = Parallelotope(A, b, [x,y])
    min_point, max_point = p.getMinPoint(), p.getMaxPoint()

    assert np.equal(max_point, [-2,-2]).all()
    assert np.equal(min_point, [-4,-4]).all()
