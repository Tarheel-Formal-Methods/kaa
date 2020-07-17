import numpy as np
import sympy as sp

from operator import add
from functools import reduce

from kaa.bernstein import BernsteinBaseConverter
from kaa.parallelotope import Parallelotope
from kaa.lputil import minLinProg, maxLinProg

from kaa.timer import Timer

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

        self.sys_dim = len(vars)
        self.num_direct = len(self.L)

    """
    Returns linear constraints representing the polytope defined by bundle.

    @returns linear constraints and their offsets.
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
    i.e each hyperplane defined by the direction vector is re-fitted to be tangent to the polytope.
    @returns canonized Bundle object
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
    @params temp_ind: index of row corresponding to desired parallelotope.i
    @returns Parallelotope object described by T[temp_ind]
    """
    def getParallelotope(self, temp_ind):

        A = np.empty([2*self.sys_dim,self.sys_dim])
        b = np.empty(2*self.sys_dim)

        'Fetch linear constraints defining the parallelotope.'
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
    Transforms the bundle according to the dynamics governing the system. (dictated by self.f)

    @params bund: Bundle object to be transformed under dynamics.
    @returns canonized transformed bundle.
    """
    def transform(self, bund):

        new_offu = np.full(bund.num_direct, np.inf)
        new_offl = np.full(bund.num_direct, np.inf)

        for row_ind, row in enumerate(bund.T):
            
            'Find the generator of the parallelotope.'
            p = bund.getParallelotope(row_ind)
            genFun = p.getGeneratorRep()

            'Create subsitutions tuples.'
            var_sub = []
            for var_ind, var in enumerate(bund.vars):
                var_sub.append((var, genFun[var_ind]))

            for column in row.astype(int):
                curr_L = bund.L[column]

                'Perform functional composition with exact transformation from unitbox to parallelotope.'
                Timer.start('Functional Composition')
                fog = [ f.subs(var_sub) for f in self.f ]
                Timer.stop('Functional Composition')

                bound_polyu = [ curr_L[func_ind] * func for func_ind, func in enumerate(fog) ]
                bound_polyu = reduce(add, bound_polyu)

                'Calculate min/max Bernstein coefficients.'
                Timer.start('Bernstein Computation')
                base_convertu = BernsteinBaseConverter(bound_polyu, bund.vars)
                max_bern_coeffu, min_bern_coeffu = base_convertu.computeBernCoeff()
                Timer.stop('Bernstein Computation')

                new_offu[column] = min(max_bern_coeffu, new_offu[column])
                new_offl[column] = min(-1 * min_bern_coeffu, new_offl[column])

        #print("New Offu: {}   NewOffl: {}".format(new_offu,new_offl))

        trans_bund = Bundle(bund.T, bund.L, new_offu, new_offl, bund.vars)
        canon_bund = trans_bund.canonize()
        return canon_bund
