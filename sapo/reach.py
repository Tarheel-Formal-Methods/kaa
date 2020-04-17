import sapo.bundle

class ReachSet:

    def __init__(self, model):
        self.model = model

    #Compute reachable set for the alloted number of time steps.
    def computeReachSet(self, time_steps):

        initial_set = self.model.bundle
        flowpipe = [initial_set]

        for ind in range(time_steps):
            trans_bund = flowpipe[ind].transform()
            flowpipe.append(trans_bund)
            print(trans_bund)

        return flowpipe
