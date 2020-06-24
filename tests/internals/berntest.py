import sympy as sp
from kaa.bernstein import BernsteinBaseConverter

def test_univar_bern_1():
    #Define symbols
    x = sp.Symbol("x")
    p = x
    b = BernsteinBaseConverter(p,[x])

    assert b.computeBernCoeff() == (1,0)

def test_univar_bern_2():
    x = sp.Symbol("x")
    p = x**2
    b = BernsteinBaseConverter(p,[x])

    assert b.computeBernCoeff() == (1,0)

def test_univar_bern_3():
    x = sp.Symbol("x")
    p = x + x**2
    b = BernsteinBaseConverter(p,[x])

    assert b.computeBernCoeff() == (2, 0)

def test_univar_bern_4():
    x = sp.Symbol("x")
    p = 3*x + 2*x**2 + x**3
    b = BernsteinBaseConverter(p,[x])

    assert b.computeBernCoeff() == (6, 0)

def test_multivar_bern_1():
    x, y = sp.Symbol("x"), sp.Symbol("y")
    p = x*y
    b = BernsteinBaseConverter(p,[x,y])

    assert b.computeBernCoeff() == (1, 0)

def test_multivar_bern_2():
    x, y = sp.Symbol("x"), sp.Symbol("y")
    p = x*y + x**2*y**2 + x**3
    b = BernsteinBaseConverter(p,[x,y])

    assert b.computeBernCoeff() == (1, 0)
