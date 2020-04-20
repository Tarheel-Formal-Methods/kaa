from sapo.bundle import Bundle, BundleTransformer
from sapo.flowpipe import FlowPipe

class ReachSet:

    def __init__(self, model):
        self.model = model

    #Compute reachable set for the alloted number of time steps.
    def computeReachSet(self, time_steps):

        initial_set = self.model.bundle
        trans = BundleTransformer(self.model.f)
        flowpipe = [initial_set]

        for ind in range(time_steps):
            trans_bund = trans.transform(flowpipe[ind])
            flowpipe.append(trans_bund)
            #print(trans_bund)

        return FlowPipe(flowpipe, self.model.vars)
