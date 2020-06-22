import numpy as np
import sympy as sp

from operator import add
from functools import reduce

from sapo.bernstein import BernsteinBaseConverter
from sapo.parallelotope import Parallelotope
from sapo.lputil import minLinProg, maxLinProg

import sapo.benchmark as Benchmark
from sapo.benchmark import Label

import sapo.log as Log
from sapo.log import Debug

class Bundle:

    def __init__(self, T, L, offu, offl, vars):

        if np.size(L,0) != np.size(offu,0):
            print("Directions matrix L and upper offsets must have matching dimensions")
            exit()

        if np.size(L,0) != np.size(offl,0):
            print("Directions matrix L and lower offsets must have matching dimensions")
            exit()

        if np.size(T,1) != np.size(L,1):
            print("Template matrix T must have the same dimensions as Directions matrix L")
            exit()

        self.T = T
        self.L = L
        self.offu = offu
        self.offl = offl
        self.vars = vars

        'Bundle statistics'
        self.sys_dim = len(vars)
        self.num_direct = len(self.L)

    """
    Returns linear constraints representing the polytope defined by their intersection.
    """
    def getIntersect(self):
        A = np.empty([2*self.num_direct,self.sys_dim])
        b = np.empty(2*self.num_direct) #row vector

        for ind in range(self.num_direct):
            A[ind] = self.L[ind]
            A[ind + self.num_direct] = np.negative(self.L[ind])
            b[ind] = self.offu[ind]
            b[ind + self.num_direct] = self.offl[ind]

        return (A,b)

    """
    Returns the bundle in its canonical form.
    """
    def canonize(self):

        A, b = self.getIntersect()

        canon_offu = np.empty(self.num_direct)
        canon_offl = np.empty(self.num_direct)

        for row_ind, row in enumerate(self.L):

            canon_offu[row_ind] = (maxLinProg(row, A, b)).fun
            canon_offl[row_ind] = (maxLinProg(np.negative(row), A, b)).fun

        return Bundle(self.T, self.L, canon_offu, canon_offl, self.vars)

    """
    Returns the Parallelotope object defined by a row in the template matrix.
    @params temp_ind: index of row corresponding to desired parallelotope.
    """
    def getParallelotope(self, temp_ind):

        A = np.empty([2*self.sys_dim,self.sys_dim])
        b = np.empty(2*self.sys_dim) #row vector

        for fac_ind, facet in enumerate(self.T[temp_ind].astype(int)):
            A[fac_ind] = self.L[facet]
            A[fac_ind + self.sys_dim] = np.negative(self.L[facet])
            b[fac_ind] = self.offu[facet]
            b[fac_ind + self.sys_dim] = self.offl[facet]

        return Parallelotope(A, b, self.vars)

    def __str__(self):
        return ''.join([str(self.offu), '  ', str(self.offl)])


class BundleTransformer:

    def __init__(self, f):
        self.f = f

    """
    Transforms the bundle according to the dyanmics governing the system. (dictated by self.f)
    @params bund: Bundle object to be transformed under dynamics.
    """
    def transform(self, bund):

        new_offu = np.full(bund.num_direct, np.inf)
        new_offl = np.full(bund.num_direct, np.inf)

        #Print bernstein computations. Separate into subroutine taking paralleltope bundle (L,offu,offl) and c^T \cdot f as parameters and outputs max/min.
        for row_ind, curr_L in enumerate(bund.L):

            max_val, min_val = self.findExtrema(bund, curr_L)
            #print("For: {0}  Bernstein Comp :Max:{1} ,  Min: {2}".format(row_ind, max_val, min_val))
            new_offu[row_ind] = min(max_val, new_offu[row_ind])
            new_offl[row_ind] = min(-1 * min_val, new_offl[row_ind])

        #Log.write_log(p_new_offu, p_new_offl, Debug.GLOBAL_BOUND)
        print("\n New Offu: {0}   New Offp = {1}\n".format(new_offu[2], new_offl[2]))
        trans_bund = Bundle(bund.T, bund.L, new_offu, new_offl, bund.vars)
        canon_bund = trans_bund.canonize()
        return canon_bund

    """
    Returns extrema of c^T \cdot f over the parallelotope bundle, P.
    @params bund: parallelotope bundle
            c: the column vector of coefficients for c^T \cdot f
    """
    def findExtrema(self, bund, c):

        'Find intersecting unitbox'
        min_coord = [ -1 * np.inf for _ in range(bund.sys_dim) ]
        max_coord = [ np.inf for _ in range(bund.sys_dim) ]


        for row_ind, row in enumerate(bund.T):
            'Calcuate minimum and maximum points of p_i'
            p = bund.getParallelotope(row_ind)
            p_min_coord = p.getMinPoint()
            p_max_coord = p.getMaxPoint()

            max_coord = [ min(p_max_coord[i], max_coord[i]) for i in range(bund.sys_dim) ]
            min_coord = [ max(p_min_coord[i], min_coord[i]) for i in range(bund.sys_dim) ]

        #Log.write_log(row_ind, p_min_coord, p_max_coord, Debug.MINMAX)
        #print("Parallelotope:  Max: {0} ,  Min: {1}".format(p_max_coord[2], p_min_coord[2]))

        'Calculate substitutions required to map unitbox over our parallelotope'
        var_sub = []
        for var_ind, var in enumerate(bund.vars):
            var_min = min_coord[var_ind] #linprog opt results
            var_max = max_coord[var_ind]

            #print((var_max,var_min))
            transf_expr = (var_max - var_min) * var + var_min
            var_sub.append((var, transf_expr))

        'compute polynomial \Lambda_i \cdot (f(v(x)))'
        bound_polyu = [ c[func_ind] * func for func_ind, func in enumerate(self.f) ]

        bound_polyu = reduce(add, bound_polyu)
        transf_bound_polyu = bound_polyu.subs(var_sub)
        #print(transf_bound_polyu)

        'Calculate min/max Bernstein coefficients'
        base_convertu = BernsteinBaseConverter(transf_bound_polyu, bund.vars)

        bern_timer = Benchmark.assign_timer(Label.BERN)
        bern_timer.start()
        max_bern_coeffu, min_bern_coeffu = base_convertu.computeBernCoeff()
        bern_timer.end()

        return max_bern_coeffu, min_bern_coeffu
