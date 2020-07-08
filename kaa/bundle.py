import numpy as np
import sympy as sp

from operator import add
from functools import reduce

from kaa.bernstein import BernsteinBaseConverter
from kaa.parallelotope import Parallelotope
from kaa.lputil import minLinProg, maxLinProg

import kaa.benchmark as Benchmark
from kaa.benchmark import Label

import kaa.log as Log
from kaa.log import Debug

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

        self.T = T # Templates
        self.L = L # Directions
        self.offu = offu
        self.offl = offl
        self.vars = vars

        self.sys_dim = len(vars)
        self.num_direct = len(self.L)

    """
    Returns linear constraints representing the polytope defined by bundle.
    """
    def getIntersect(self):
        A = np.empty([2*self.num_direct,self.sys_dim])
        b = np.empty(2*self.num_direct) #row vector

        for ind in range(self.num_direct):
            A[ind] = self.L[ind]
            A[ind + self.num_direct] = np.negative(self.L[ind])
            b[ind] = self.offu[ind]
            b[ind + self.num_direct] = self.offl[ind]

        return A, b

    """
    Returns the bundle with tightest offsets for each direction vector in self.L
    i.e each hyperplane defined by the direction vector is tangent to the polytope.
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

        # Fetch linear constraints defining the parallelotope.
        for fac_ind, facet in enumerate(self.T[temp_ind].astype(int)):
            A[fac_ind] = self.L[facet]
            A[fac_ind + self.sys_dim] = np.negative(self.L[facet])
            b[fac_ind] = self.offu[facet]
            b[fac_ind + self.sys_dim] = self.offl[facet]

        return Parallelotope(A, b, self.vars)

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

        for row_ind, row in enumerate(bund.T):
            #Calcuate minimum and maximum points of p_i
            p = bund.getParallelotope(row_ind)
            genFun = p.getGeneratorRep()
            
            #Calculate transformation subsitutions
            var_sub = []
            for var_ind, var in enumerate(bund.vars):
                var_sub.append((var, genFun[var_ind]))

            for column in row.astype(int):
                curr_L = bund.L[column] #row in L

                #compute polynomial \Lambda_i \cdot (f(v(x)))

                bound_polyu = [ curr_L[func_ind] * func for func_ind, func in enumerate(self.f) ]

                bound_polyu = reduce(add, bound_polyu) #transform to range over unit box
                transf_bound_polyu = bound_polyu.subs(var_sub)

                #Calculate min/max Bernstein coefficients
                base_convertu = BernsteinBaseConverter(transf_bound_polyu, bund.vars)
                max_bern_coeffu, min_bern_coeffu = base_convertu.computeBernCoeff()

                #Log.write_log(max_bern_coeffu, min_bern_coeffu, row, Debug.LOCAL_BOUND)

                new_offu[column] = min(max_bern_coeffu, new_offu[column])
                new_offl[column] = min(-1 * min_bern_coeffu, new_offl[column])

        #Log.write_log(p_new_offu, p_new_offl, Debug.GLOBAL_BOUND)

        trans_bund = Bundle(bund.T, bund.L, new_offu, new_offl, bund.vars)
        canon_bund = trans_bund.canonize()
        return canon_bund

    """
    Returns extrema of c^T \cdot f over the parallelotope bundle, P.
    @params bund: parallelotope bundle
            c: the column vector of coefficients for c^T \cdot f
    """
    def findExtrema(self, bund, c):

        'Find the bounding box over the intersections of the parallelotopes in bundle.'
        min_coord = [ -1 * np.inf for _ in range(bund.sys_dim) ]
        max_coord = [ np.inf for _ in range(bund.sys_dim) ]

        for row_ind, row in enumerate(bund.T):

            'Calcuate minimum and maximum points of p_i'
            p = bund.getParallelotope(row_ind)
            p_min_coord = p.getMinPoint()
            p_max_coord = p.getMaxPoint()

            max_coord = [ min(p_max_coord[i], max_coord[i]) for i in range(bund.sys_dim) ]
            min_coord = [ max(p_min_coord[i], min_coord[i]) for i in range(bund.sys_dim) ]

        'Calculate substitutions required to map unitbox over our parallelotope'
        var_sub = []
        for var_ind, var in enumerate(bund.vars):
            var_min = min_coord[var_ind] # LP results
            var_max = max_coord[var_ind]

            transf_expr = (var_max - var_min) * var + var_min
            var_sub.append((var, transf_expr))

        'compute polynomial \Lambda_i \cdot (f(v(x)))'
        bound_polyu = [ c[func_ind] * func for func_ind, func in enumerate(self.f) ]

        bound_polyu = reduce(add, bound_polyu)
        transf_bound_polyu = bound_polyu.subs(var_sub)

        'Calculate min/max Bernstein coefficients over the calculated polynomial'
        base_convertu = BernsteinBaseConverter(transf_bound_polyu, bund.vars)
        max_bern_coeffu, min_bern_coeffu = base_convertu.computeBernCoeff()

        return max_bern_coeffu, min_bern_coeffu
