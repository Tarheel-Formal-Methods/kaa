import sympy as sp

from functools import reduce
from math import factorial
from operator import mul,add

class BernsteinBaseConverter:

    def __init__(self, poly, vars):
        self.poly = sp.Poly(poly, vars)
        self.coeffs = self.poly.coeffs()
        self.vars = vars
        self.var_num = len(vars)
        self.degree = self._getDegree()
        self.monom_list = self._computeLowerDeg(self.degree)

    #Computes the maximum and minimum Bernstein coefficients for self.poly
    def computeBernCoeff(self):

        bern_coeff = []

        for monom in self.monom_list:
            bern_coeff.append(self._computeIthBernCoeff(monom))

        return max(bern_coeff), min(bern_coeff)

    #Compute the ith Bernstein coefficient.
    def _computeIthBernCoeff(self, i):

        lower_degs = self._computeLowerDeg(i) #redundent operation

        bern_sum_list = []
        for j in lower_degs:
            coeff_mul_list = []

            for ind, _ in enumerate(self.degree):
                j_coeff = self._choose(i[ind], j[ind]) / self._choose(self.degree[ind], j[ind])
                coeff_mul_list.append(j_coeff)

            bern_coef = reduce(mul,coeff_mul_list)
            poly_coef = self.poly.coeff_monomial(self._getMonom(j))
            #print(i,poly_coef)
            bern_sum_list.append(bern_coef * poly_coef)

        #print(bern_sum_list)
        return reduce(add, bern_sum_list)

    #Compute the ith Bernstein polynomial with degree tuple i
    def _computeIthBernBasis(self, i):

        bern_poly_list = []
        for var, var_index in enumerate(self.vars):
            curr_deg = self.degree[var_index]
            curr_i = i[var_index]

            ith_coeff = self._choose(curr_deg, curr_i)
            ith_bern_expr = ith_coeff * var**curr_i + (1 - var)**(curr_deg - curr_i)
            bern_poly_list.append(ith_bern_expr)

        return functools.reduce(mul, bern_poly_list)

    def _computeLowerDeg(self, i):

        assert len(i) == len(self.degree)

        degree_list = []
        accum_list = [0 for _ in i]

        def recurse(curr_degree, degree_pos):
            if degree_pos == len(i) - 1:
                degree_list.append(curr_degree.copy())
                if curr_degree[degree_pos] < i[degree_pos]:
                    curr_degree[degree_pos] += 1
                    recurse(curr_degree, degree_pos)
            else:
                for _ in range(i[degree_pos]+1):
                    recurse(curr_degree, degree_pos + 1)
                    curr_degree[degree_pos] += 1
                    curr_degree[degree_pos + 1] = 0

        recurse(accum_list, 0)
        return degree_list

    def _getDegree(self):

        monom_tups = self.poly.monoms()
        degree = []

        for var_index, _ in enumerate(self.vars):
             var_deg = max([ monom[var_index] for monom in monom_tups])
             degree.append(var_deg)

        return degree

    #Get monomial indexed by degree tuple
    def _getMonom(self, j):

        var_monom = [self.vars[i]**j[i] for i in range(self.var_num)]
        expr = reduce(mul, var_monom)

        monomial = sp.Poly(expr, self.vars)

        return monomial.as_expr()

    def _choose(self, n, r):
        return factorial(n) / (factorial(r)*factorial(n-r))
