import sympy as sp

from bernstein import *


def main():
    x, y = sp.Symbol("x"), sp.Symbol("y")
    p = 2*x*y + 5*y**2 + x**2

    b = BernsteinBaseConverter([x,y], p)
    print(b.computeBernCoeff())

if __name__ == "__main__":
    main()
