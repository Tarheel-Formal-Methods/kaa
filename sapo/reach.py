import sapo.benchmark as Benchmark
import sapo.log as Log

from sapo.benchmark import Label
from sapo.log import Debug
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

            transf_timer = Benchmark.assign_timer(Label.TRANSF)

            #Log.write_log(ind, Debug.STEP)
            transf_timer.start()
            trans_bund = trans.transform(flowpipe[ind])
            transf_timer.end()

            print("Computed Step {0} -- Time Elapsed: {1} sec".format(ind, transf_timer.duration))
            flowpipe.append(trans_bund)

        return FlowPipe(flowpipe, self.model.vars)
