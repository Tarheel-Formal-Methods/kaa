import numpy as np
import sympy as sp
from scipy.optimize import linprog

from operator import add
from functools import reduce

from sapo.bernstein import BernsteinBaseConverter
from sapo.parallelotope import Parallelotope

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
        self.vars = vars #List of variables used in bundle computation

        #bundle stats
        self.sys_dim = len(vars)
        self.num_direct = len(self.L)

    #Get the intersection of the polytope
    def getIntersect(self):
        A = np.empty([2*self.num_direct,self.sys_dim])
        b = np.empty(2*self.num_direct) #row vector

        for ind in range(self.num_direct):
            A[ind] = self.L[ind]
            A[ind + self.num_direct] = -1 * self.L[ind]
            b[ind] = self.offu[ind]
            b[ind + self.num_direct] = self.offl[ind]

        return (A,b)

    #Canonize the bundle.
    def canonize(self):

        A, b = self.getIntersect()

        canon_offu = np.empty(self.num_direct)
        canon_offl = np.empty(self.num_direct)

        for row_ind, row in enumerate(self.L):

            canon_offu[row_ind] = (linprog(row, A, b)).fun
            canon_offl[row_ind] = (linprog(np.negative(row), A, b)).fun

        return Bundle(self.T, self.L, canon_offu, canon_offl, self.vars)

    #get Parallelotope object associated to template_ind
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

    # Transform the bundle according to Polynomial DDS
    def transform(self, bund):

        p_new_offu = np.full(bund.num_direct, np.inf)
        p_new_offl = np.full(bund.num_direct, np.inf)

        #get parallelotope P_i
        for row_ind, row in enumerate(bund.T):

            #Calcuate minimum and maximum points of p_i
            p = bund.getParallelotope(row_ind)
            p_min_coord = p.getMinPoint()
            p_max_coord = p.getMaxPoint()

            #Calculate transformation subsitutions
            var_sub = []
            for var_ind, var in enumerate(bund.vars):
                p_min = p_min_coord[var_ind] #linprog opt results
                p_max = p_max_coord[var_ind]

                transf_expr = (p_max - p_min) * var + p_min
                var_sub.append((var, transf_expr))

            for column in row.astype(int):
                #get facet
                curr_L = bund.L[column] #row in L

                #compute polynomial \Lambda_i \cdot (f(v(x)))

                bound_polyu = [ curr_L[func_ind] * func for func_ind, func in enumerate(self.f) ]
                bound_polyl = [ -1 * curr_L[func_ind] * func for func_ind, func in enumerate(self.f) ]

                bound_polyu = reduce(add, bound_polyu)
                bound_polyl = reduce(add, bound_polyl)

                #transform to range over unit box
                transf_bound_polyu = bound_polyu.subs(var_sub)
                transf_bound_polyl = bound_polyl.subs(var_sub)
                #print(transf_bound_polyu, transf_bound_polyl)

                #Calculate min/max Bernstein coefficients
                base_convertu = BernsteinBaseConverter(transf_bound_polyu, bund.vars)
                max_bern_coeffu, min_bern_coeffu = base_convertu.computeBernCoeff()

                base_convertl = BernsteinBaseConverter(transf_bound_polyl, bund.vars)
                max_bern_coeffl, min_bern_coeffl = base_convertl.computeBernCoeff()

                #print(''.join(['Max:', str(max_bern_coeffu),' Min: ', str(max_bern_coeffl), 'for P: ', str(row)]))
                p_new_offu[column] = min(max_bern_coeffu, p_new_offu[column])
                p_new_offl[column] = min(max_bern_coeffl, p_new_offl[column])

        trans_bund = Bundle(bund.T, bund.L, p_new_offu, p_new_offl, bund.vars) #Major issues could arise with unused direcitions
        canon_bund = trans_bund.canonize()
        return canon_bund
