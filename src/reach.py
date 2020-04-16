from . import bundle

class ReachSet:


def __init__(self, model):
    self.model = model

#Compute reachable set for the alloted number of time steps.
def computeReachSet(self, time_steps):

    initial_set = model.bundle
    flowpipe = [].append(initial_set)

    for ind in range(time_steps):

        trans_bund = flowpipe[ind].transform()
        flowpipe.append(trans_bund)

    return flowpipe
