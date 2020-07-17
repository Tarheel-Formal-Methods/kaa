from kaa.timer import Timer
import kaa.log as Log

from kaa.log import Debug
from kaa.bundle import Bundle, BundleTransformer
from kaa.flowpipe import FlowPipe

'''
Object handling all reachable flowpipe computations.
'''
class ReachSet:

    def __init__(self, model):
        self.model = model

    """
    Compute reachable set for the alloted number of time steps.
    @params time_steps: number of time steps to carry out the reachable set computation.
    @returns FlowPipe object with computed flowpipe
    """
    def computeReachSet(self, time_steps):

        initial_set = self.model.bund
        transformer = BundleTransformer(self.model.f)
        flowpipe = [initial_set]


        for ind in range(time_steps):


            Timer.start('Reachable Set Computation')
            trans_bund = transformer.transform(flowpipe[ind])
            reach_time = Timer.stop('Reachable Set Computation')

            print("Computed Step {0} -- Time Elapsed: {1} sec".format(ind, reach_time))
            flowpipe.append(trans_bund)

        return FlowPipe(flowpipe, self.model.vars)
