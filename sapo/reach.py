import time

from sapo.bundle import Bundle, BundleTransformer
from sapo.flowpipe import FlowPipe

"""
Object handling all reachable flowpipe computations.
TODO
    - Proper debugging framework.
"""
class ReachSet:

    def __init__(self, model):
        self.model = model

    #Compute reachable set for the alloted number of time steps.
    def computeReachSet(self, time_steps):

        initial_set = self.model.bundle
        trans = BundleTransformer(self.model.f)
        flowpipe = [initial_set]

        for ind in range(time_steps):

            start = time.time()
            trans_bund = trans.transform(flowpipe[ind])
            end = time.time()

            print("Computed Step {0} -- Time Elapsed: {1} sec".format(ind, end - start))
            flowpipe.append(trans_bund)

        return FlowPipe(flowpipe, self.model.vars)
