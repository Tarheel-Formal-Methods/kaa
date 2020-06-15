import sympy as sp

from functools import reduce
from math import factorial
from operator import mul,add
from itertools import product

class BernsteinBaseConverter:

    def __init__(self, poly, vars):
        self.poly = sp.Poly(poly, vars)
        self.coeffs = self.poly.coeffs()
        self.vars = vars
        self.var_num = len(vars)
        self.degree = self._getDegree()
        self.monom_list = self._computeLowerDeg(self.degree)

    """
    Computes the maximum and minimum Bernstein coefficients for self.poly
    Notes: This can be most likely improved through an improved method through matrices.
    """
    def computeBernCoeff(self):

        bern_coeff = []

        for monom in self.monom_list:
            bern_coeff.append(self._computeIthBernCoeff(monom))

        return max(bern_coeff), min(bern_coeff)

    """
    Compute the ith Bernstein coefficient.
    """
    def _computeIthBernCoeff(self, i):

        lower_degs = self._computeLowerDeg(i)

        bern_sum_list = []
        for j in lower_degs:
            coeff_mul_list = []

            for ind, _ in enumerate(self.degree):
                j_coeff = self._choose(i[ind], j[ind]) / self._choose(self.degree[ind], j[ind])
                coeff_mul_list.append(j_coeff)

            bern_coef = reduce(mul,coeff_mul_list)
            poly_coef = self.poly.coeff_monomial(self._getMonom(j))
            bern_sum_list.append(bern_coef * poly_coef)

        return reduce(add, bern_sum_list)

    def _computeLowerDeg(self, i):

        assert len(i) == len(self.degree)

        iterators = [range(idx+1) for idx in i]
        return list(product(*iterators))

    def _getDegree(self):

        monom_tups = self.poly.monoms()
        degree = []
        
        for var_index, _ in enumerate(self.vars):
             var_deg = max([ monom[var_index] for monom in monom_tups])
             degree.append(var_deg)

        return degree

    def _getMonom(self, j):

        var_monom = [self.vars[i]**j[i] for i in range(self.var_num)]
        expr = reduce(mul, var_monom)

        monomial = sp.Poly(expr, self.vars)
        return monomial.as_expr()

    def _choose(self, n, r):
        return factorial(n) // (factorial(r)*factorial(n-r))
