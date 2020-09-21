import numpy as np
import multiprocessing as mp

from kaa.lputil import minLinProg, maxLinProg
from kaa.timer import Timer
from kaa.settings import KaaSettings

"""
Object encapsulating routines calculating properties of parallelotopes.
"""
class Parallelotope:

    def __init__(self, A, b, vars):
        self.vars = vars
        self.dim = len(vars)

        self.A = A[:self.dim]
        self.b = b

    """
    Return list of functions transforming the n-unit-box over the parallelotope.
    @returns list of transfomation from unitbox over the parallelotope.
    """
    def getGeneratorRep(self):

        Timer.start('Generator Procedure')
        base_vertex = self._computeBaseVertex()
        gen_list = self._computeGenerators(base_vertex)

        'Create list representing the linear transformation q + \sum_{j} a_j* g_j'
        expr_list = base_vertex
        for var_ind, var in enumerate(self.vars):
            for i in range(self.dim):
                expr_list[i] += gen_list[var_ind][i] * var
        Timer.stop('Generator Procedure')

        return expr_list
    
    """
    Calculate generators as substraction: vertices - base_vertex.
    We calculate the vertices by solving the following linear system for each vertex i:


    Ax = [b_1, ... , -b_{i+n}, ... , b_n]^T

    Note that this simply finds the vertex to calculate the generator vectors. The generators will be the vectors g_j = v_j - q
    where v_j is the jth vertex calculated in the manner above.

    The parallelotope will be exprssed as sum of the base vertex and the convex combination of the generators.

    p(a_1, ... ,a_n) =  q + \sum_{j} a_j * g_j

    where q is the base vertex and the g_j are the generators. a_j will be in the unitbox [0,1]

    @params base_vertex: base vertex q
    @returns generator vectors g_j
    """
    def _computeGenerators(self, base_vertex):

        u_b = self.b[:self.dim]
        coeff_mat = self.A

        'Hacky way to toggle parallelism for experiments'
        if KaaSettings.use_parallel:
            p = mp.Pool(processes=4)
            vertices = p.starmap(self._gen_worker, [ (i, u_b, coeff_mat) for i in range(self.dim) ])
            p.close()
            p.join()
        else:
            vertices = []
            for i in range(self.dim):
               vertices.append(self._gen_worker(i, u_b, coeff_mat))

        return [ [x-y for x,y in zip(vertices[i], base_vertex)] for i in range(self.dim) ]

    """
    Worker process for calculating vertices of higher-dimensional parallelotopes.
    Only called by Pool.starmap
    @params i - vertex index
            u_b, coef_mat - shared reference to upper offsets and directions matrix.
    @returns coordinates of vertex
    """
    def _gen_worker(self, i, u_b, coeff_mat):

        negated_bi = np.copy(u_b)
        negated_bi[i] = -self.b[i + self.dim]
        sol_set_i = np.linalg.solve(self.A, negated_bi)

        return sol_set_i


    """
    Calculate the base vertex of the parallelotope (variable q)
    We calculate the vertices by solving a linear system of the following form:

    Ax = u_b

    where u_b are the offsets for the upper facets of parallelotope (first half of self.b).

    @returns base-vertex in list
    """
    def _computeBaseVertex(self):

        u_b = self.b[:self.dim]

        sol_set = np.linalg.solve(self.A, u_b)
        return list(sol_set)
