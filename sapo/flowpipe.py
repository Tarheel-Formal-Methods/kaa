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
        self.dim_sys = len(flowpipe)

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

    def plot2DPhase(self, var_xind, var_yind):

        #Define the following projected normal vectors
        x,y = flowpipe.vars[var_xind], flowpipe.vars[var_yind]
        norm_vecs = [[0 for _ in range(dim_sys)] for _ in range(4)]

        norm_vecs[0][x] = 1
        norm_vecs[1][y] = 1
        norm_vecs[2][x] = -1
        norm_vecs[3][y] = -1

        norm_A = np.stack(norm_vecs, axis=0)

        fig, ax = plt.subplots(1,1)

        #Effort into renaming some of these vars
        offsets = []
        c = [0, 0, -1]
        for bund in self.flowpipe:
            bund_A, bund_b = bund.getIntersect()

            bund_off = np.empty(4)
            for i in range(4):
                bund_off[i] = -1 * linprog(np.negative(norm_vecs[i]), bund_A, bund_b, bounds=(None,None)).fun

            phase_intersect = np.stack(norm_A,bund_off)
            #compute center of intersection

            center_A = np.empty([4,3])
            for k in range(4):
                center_A[k] = np.concatenate(norm_vecs[k], [np.pdist(norm_vecs[k], 'euclidean')])

            center_pt = linprog(c, center_A, bund_off, bounds=(None,None)).x

            hs = HalfspaceIntersection(phase_intersect, center_pt)
            x, y = zip(*hs.intersections)
            ax.fill(x, y)
        fig,show()
