import numpy as np
from sympy import Matrix, linsolve, EmptySet

from kaa.lputil import minLinProg, maxLinProg

"""
Object encapsulating routines calculating properties of parallelotopes
"""
class Parallelotope:

    def __init__(self, A, b, vars):
        self.vars = vars
        self.dim = len(vars)

        self.A = A[:self.dim] #Tailored to find base vertex and generators.
        self.b = b
        self.d = b[:self.dim]

    '''
    Return list of functions transforming the n-unit-box over the parallelotope.
    '''

    def getGeneratorRep(self):

        base_vertex = self._computeBaseVertex()
        gen_list = self._computeGenerators(base_vertex)
        #print("Gen List:", gen_list)

        expr_list = base_vertex
        for var_ind, var in enumerate(self.vars):
            for i in range(self.dim):
                expr_list[i] = expr_list[i] + gen_list[var_ind][i] * var

        #print("Expr List:" , expr_list)
        return expr_list
    
    '''
    Calculate generators as substraction: vertices - base_vertex
    @params base_vertex: base vertex
    '''
    def _computeGenerators(self, base_vertex):

        vertices = []
        coeff_mat = self._convertMatFormat(self.A)
        for i in range(self.dim):
            negated_bi = np.copy(self.d)
            negated_bi[i] = -self.b[i + self.dim]
            negated_bi = self._convertMatFormat(negated_bi)
            #print("Negated bi:", negated_bi)
            sol_set_i = linsolve((coeff_mat, negated_bi), self.vars)
            vertex_i = self._convertSolSetToList(sol_set_i)
            vertices.append(vertex_i)

        #print("Vertices:", vertices)
        return [ [x-y for x,y in zip(vertices[i], base_vertex)] for i in range(self.dim) ]
        

    '''
    Calculate the base vertex of the parallelotope (variable q)
    '''
    def _computeBaseVertex(self):
        coeff_mat = self._convertMatFormat(self.A)
        offset_mat = self._convertMatFormat(self.d)

        #print('A: {}   b: {}'.format(coeff_mat,offset_mat))
        sol_set = linsolve((coeff_mat, offset_mat), self.vars)

        return self._convertSolSetToList(sol_set)

    '''
    Convert numpy matrix into sympy matrix
    @params mat: numpy matrix
    '''
    def _convertMatFormat(self, mat):
        rows = mat.tolist()
        return Matrix(rows)
    '''
    Takes solution set returned by sympy and converts into list
    @params fin_set: FiniteSet
    '''
    def _convertSolSetToList(self, fin_set):

        assert fin_set is not EmptySet

        return list(fin_set.args[0])
