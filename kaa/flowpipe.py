import matplotlib.pyplot as plt
import numpy as np

from scipy.spatial import HalfspaceIntersection

import kaa.log as Log
from kaa.log import Debug
from kaa.lputil import minLinProg, maxLinProg
from kaa.timer import Timer

"""
Object encapsulating flowpipe data.
"""
class FlowPipe:

    def __init__(self, flowpipe, vars):

        self.flowpipe = flowpipe
        self.vars = vars

    def __len__(self):
        return len(self.flowpipe)

    def __iter__(self):
        return iter(self.flowpipe)

"""
Plotting object handling all matplotlib manipulations and plot-formatting.
"""
class FlowPipePlotter:

    def __init__(self, flowpipe):
        self.flowpipe = flowpipe
        self.dim_sys = len(flowpipe.vars)
        self.vars = self.flowpipe.vars

    """
    Plots projection of reachable set against time t.
    @params (vars_tup): tuple containing indices of variables to be plotted.

    Only plots to an on-screen figure in its current state.
    """

    def plot2DProj(self, *vars_tup):
        num_var = len(vars_tup)
        pipe_len = len(self.flowpipe)

        fig, ax = plt.subplots(1, num_var)
        ax = [ax] if num_var == 1 else ax

        'Time step vector for plotting'
        t = np.arange(0, pipe_len, 1)

        for ax_ind, var_ind in enumerate(vars_tup):

            Timer.start('Proj')
            curr_var = self.vars[var_ind]

            'Vector of minimum and maximum points of the polytope represented by parallelotope bundle.'
            y_min, y_max = np.empty(pipe_len), np.empty(pipe_len)

            'Initialize objective function'
            y_obj = [0 for _ in self.vars]
            y_obj[var_ind] = 1

            'Calculate the minimum and maximum points through LPs for every iteration of the bundle.'
            for bund_ind, bund in enumerate(self.flowpipe):

                bund_A, bund_b = bund.getIntersect()

                y_min[bund_ind] = minLinProg(y_obj, bund_A, bund_b).fun
                y_max[bund_ind] = maxLinProg(y_obj, bund_A, bund_b).fun

            ax[ax_ind].fill_between(t, y_min, y_max)
            ax[ax_ind].set_xlabel("t: time steps")
            ax[ax_ind].set_ylabel("Reachable Set for {}".format(curr_var))

            plot_time = Timer.stop('Proj')
            print("Plotting projection for dimension {} done -- Time Spent: {}".format(curr_var, plot_time))

        fig.show()

    """
    Plots phase between two variables of dynamical system.
    @params x: index of variable to be plotted as x-axis of desired phase
            y: index of variable to be plotted as y-axis of desired phase
    """
    def plot2DPhase(self, x, y):


        Timer.start('Phase')

        x_var, y_var = self.vars[x], self.vars[y]

        'Define the following projected normal vectors.'
        norm_vecs = np.zeros([6,self.dim_sys])
        norm_vecs[0][x] = 1; norm_vecs[1][y] = 1;
        norm_vecs[2][x] = -1; norm_vecs[3][y] = -1;
        norm_vecs[4][x] = 1; norm_vecs[4][y] = 1;
        norm_vecs[5][x] = -1; norm_vecs[5][y] = -1;

        fig, ax = plt.subplots(1)
        comple_dim = [i for i in range(self.dim_sys) if i not in [x,y]]

        'Initialize objective function for Chebyshev intersection LP routine.'
        c = [0 for _ in range(self.dim_sys + 1)]
        c[-1] = 1

        for bund in self.flowpipe:
            bund_A, bund_b = bund.getIntersect()

            'Compute the normal vector offsets'
            bund_off = np.empty([len(norm_vecs),1])
            for i in range(len(norm_vecs)):
                bund_off[i] = minLinProg(np.negative(norm_vecs[i]), bund_A, bund_b).fun

            'Remove irrelevant dimensions. Mostly doing this to make HalfspaceIntersection happy.'
            phase_intersect = np.hstack((norm_vecs, bund_off))
            phase_intersect = np.delete(phase_intersect, comple_dim, axis=1)

            'Compute Chebyshev center of intersection.'
            row_norm = np.reshape(np.linalg.norm(norm_vecs, axis=1), (norm_vecs.shape[0],1))
            center_A = np.hstack((norm_vecs, row_norm))

            neg_bund_off = np.negative(bund_off)
            center_pt = maxLinProg(c, center_A, list(neg_bund_off.flat)).x
            center_pt = np.asarray([b for b_i, b in enumerate(center_pt) if b_i in [x, y]])

            'Run scipy.spatial.HalfspaceIntersection.'
            hs = HalfspaceIntersection(phase_intersect, center_pt)
            inter_x, inter_y = zip(*hs.intersections)
            ax.set_xlabel('{}'.format(x_var))
            ax.set_ylabel('{}'.format(y_var))
            ax.fill(inter_x, inter_y, 'b')

        fig.show()

        phase_time = Timer.stop('Phase')
        print("Plotting phase for dimensions {}, {} done -- Time Spent: {}".format(x_var, y_var, phase_time))
