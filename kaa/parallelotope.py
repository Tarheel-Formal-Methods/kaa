import numpy as np
from sympy import Matrix, linsolve, EmptySet

from kaa.lputil import minLinProg, maxLinProg
from kaa.timer import Timer

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

    Note that this simply finds the vertex to calculate the generator vectors.
    The parallelotope will be exprssed as sum of the base vertex and the convex combination of the generators.

   p(a_1, ... ,a_n) =  q + \sum_{j} a_j * g_j

    where q is the base vertex and the g_j are the generators. a_j will be in the unitbox [0,1]

    @params base_vertex: base vertex q
    @returns generator vectors g_j
    """
    def _computeGenerators(self, base_vertex):

        u_b = self.b[:self.dim]

        vertices = []
        coeff_mat = self._convertMatFormat(self.A)

        for i in range(self.dim):
            negated_bi = np.copy(u_b)
            negated_bi[i] = -self.b[i + self.dim]
            negated_bi = self._convertMatFormat(negated_bi)
            
            sol_set_i = linsolve((coeff_mat, negated_bi), self.vars)
            vertex_i = self._convertSolSetToList(sol_set_i)
            vertices.append(vertex_i)

        return [ [ x-y for x,y in zip(vertices[i], base_vertex) ] for i in range(self.dim) ]
       

    """
    Calculate the base vertex of the parallelotope (variable q)
    We calculate the vertices by solving a linear system of the following form:

    Ax = u_b

    where u_b are the offsets for the upper facets of parallelotope (first half of self.b).
    @returns base-vertex in list
    """
    def _computeBaseVertex(self):

        u_b = self.b[:self.dim]

        coeff_mat = self._convertMatFormat(self.A)
        offset_mat = self._convertMatFormat(u_b)

        sol_set = linsolve((coeff_mat, offset_mat), self.vars)
        return self._convertSolSetToList(sol_set)

    """
    Convert numpy matrix into sympy matrix
    @params mat: numpy matrix
    @returns sympy matrix counterpart
    """
    def _convertMatFormat(self, mat):
        return Matrix(mat.tolist())
    
    """
    Takes solution set returned by sympy and converts into list
    @params fin_set: FiniteSet
    @returns list of sympy solution set
    """
    def _convertSolSetToList(self, fin_set):

        assert fin_set is not EmptySet

        return list(fin_set.args[0])
