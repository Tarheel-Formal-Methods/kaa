import matplotlib.pyplot as plt
import numpy as np

from scipy.optimize import linprog

class FlowPipe:

    def __init__(self, flowpipe, vars):

        self.flowpipe = flowpipe
        self.vars = vars

    def __len__(self):
        return len(self.flowpipe)

    def __iter__(self):
        return iter(self.flowpipe)

class FlowPipePlotter:

    def __init__(self, flowpipe):
        self.flowpipe = flowpipe

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
            y_min[bund_ind] = (linprog(y_min_obj, bund_A, bund_b)).fun
            y_max[bund_ind] = -1 * (linprog(y_max_obj, bund_A, bund_b)).fun

        #print(y_min, y_max)
        ax.fill_between(t, y_min, y_max)
        ax.set_xlabel("t: time steps")
        ax.set_ylabel("var: Reachable Set of Chosen Dim")
        fig.show()
