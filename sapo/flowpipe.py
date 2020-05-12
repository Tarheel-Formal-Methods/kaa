import matplotlib.pyplot as plt
import numpy as np

from scipy.optimize import linprog
from scipy.spatial import HalfspaceIntersection

class FlowPipe:

    def __init__(self, flowpipe, vars):

        self.flowpipe = flowpipe
        self.vars = vars

    def __len__(self):
        return len(self.flowpipe)

    def __iter__(self):
        return iter(self.flowpipe)

"""
Plotting object handling all matplotlib tools and plot formatting.
TODO
    - Ability to plot several projections and phase plots against each other neatly.
"""
class FlowPipePlotter:

    def __init__(self, flowpipe):
        self.flowpipe = flowpipe
        self.dim_sys = len(flowpipe.vars)

    """
    Plots projection of reachable set against time t.
    """

    def plot2DProj(self, var_ind):
        fig, ax = plt.subplots(1,1)
        pipe_len = len(self.flowpipe)

        t = np.arange(0, pipe_len, 1)
        y_min = np.empty(pipe_len)
        y_max = np.empty(pipe_len)

        y_min_obj = [0 for _ in self.flowpipe.vars]
        y_min_obj[var_ind] = 1

        y_max_obj = [0 for _ in self.flowpipe.vars]
        y_max_obj[var_ind] = -1

        for bund_ind, bund in enumerate(self.flowpipe):
            bund_A, bund_b = bund.getIntersect()
            y_min[bund_ind] = (linprog(y_min_obj, bund_A, bund_b, bounds=(None,None))).fun
            y_max[bund_ind] = -1 * (linprog(y_max_obj, bund_A, bund_b, bounds=(None,None))).fun

        ax.fill_between(t, y_min, y_max)
        ax.set_xlabel("t: time steps")
        ax.set_ylabel("Reachable Set for {0}".format(self.flowpipe.vars[var_ind]))
        fig.show()

    def plot2DPhase(self, x, y):

        #Define the following projected normal vectors
        norm_vecs = np.zeros([8,self.dim_sys])

        norm_vecs[0][x] = 1; norm_vecs[1][y] = 1; #Testing support functions for these normals for now
        norm_vecs[2][x] = -1; norm_vecs[3][y] = -1;
        norm_vecs[4][x] = 1;  norm_vecs[4][y] = 1;
        norm_vecs[5][x] = 1;  norm_vecs[5][y] = -1;
        norm_vecs[6][x] = -1;  norm_vecs[6][y] = 1;
        norm_vecs[7][x] = -1;  norm_vecs[7][y] = -1;

        fig, ax = plt.subplots(1,1)

        #Effort into renaming some of these vars
        comple_dim = [i for i in range(self.dim_sys) if i not in [x,y]] #dimensions not of x,y

        c = [0 for _ in range(self.dim_sys + 1)]
        c[-1] = -1

        for bund in self.flowpipe:
            bund_A, bund_b = bund.getIntersect()

            bund_off = np.empty([len(norm_vecs),1])
            for i in range(len(norm_vecs)):
                bund_off[i] = linprog(np.negative(norm_vecs[i]), bund_A, bund_b, bounds=(None,None)).fun

            phase_intersect = np.hstack((norm_vecs, bund_off))  #remove irrelevant dimensions. Mostly doing this to make HalfspaceIntersection happy.
            phase_intersect = np.delete(phase_intersect, comple_dim, axis=1)

            #compute center of intersection
            row_norm = np.reshape(np.linalg.norm(norm_vecs, axis=1), (norm_vecs.shape[0],1))
            center_A = np.hstack((norm_vecs, row_norm))

            #print(np.negative(bund_off.T))
            center_pt = linprog(c, A_ub=center_A, b_ub=np.negative(bund_off.T)).x
            center_pt = np.asarray([b for b_i, b in enumerate(center_pt) if b_i in [x, y]])

            hs = HalfspaceIntersection(phase_intersect, center_pt) #Issues with facets and feasible point being coplanar.
            inter_x, inter_y = zip(*hs.intersections)
            ax.fill(inter_x, inter_y, 'b')

        fig.show()
