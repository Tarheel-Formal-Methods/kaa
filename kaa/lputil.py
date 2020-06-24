'''
Wrapper interface to GLPK
'''
import swiglpk as glpk
from itertools import product

class LPSolution:

    def __init__(self, x, fun):
        self.x = x
        self.fun = fun

minLinProg = lambda c, A, b: _linprog(c, A, b, glpk.GLP_MIN)
maxLinProg = lambda c, A, b: _linprog(c, A, b, glpk.GLP_MAX)

def _linprog(c, A, b, obj):

    lp = glpk.glp_create_prob()
    glpk.glp_set_obj_dir(lp, obj)

    params = glpk.glp_smcp()
    glpk.glp_init_smcp(params)
    params.msg_lev = glpk.GLP_MSG_OFF #Only print error messages from GLPK

    num_rows = A.shape[0]
    num_cols = A.shape[1]
    mat_size = num_rows * num_cols

    glpk.glp_add_rows(lp, num_rows)

    for row_ind in range(num_rows):
        glpk.glp_set_row_bnds(lp, row_ind+1, glpk.GLP_UP, 0.0, float(b[row_ind]))

    glpk.glp_add_cols(lp, num_cols)

    for col_ind in range(num_cols):
        glpk.glp_set_col_bnds(lp, col_ind+1, glpk.GLP_FR, 0.0, 0.0)
        glpk.glp_set_obj_coef(lp, col_ind+1, c[col_ind])

    'Swig arrays are used for feeding constraints in GLPK'

    ia, ja, ar = [],[],[]
    for i,j in product(range(num_rows),range(num_cols)):
        ia.append(i+1)
        ja.append(j+1)
        ar.append(float(A[i][j]))

    ia = glpk.as_intArray(ia)
    ja = glpk.as_intArray(ja)
    ar = glpk.as_doubleArray(ar)

    glpk.glp_load_matrix(lp, mat_size, ia, ja, ar)
    glpk.glp_simplex(lp, params)

    fun = glpk.glp_get_obj_val(lp)
    x = [i for i in map(lambda x: glpk.glp_get_col_prim(lp, x+1), range(num_cols))]

    glpk.glp_delete_prob(lp)
    glpk.glp_free_env()

    return LPSolution(x, fun)
