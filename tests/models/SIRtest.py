from sapo.reach import ReachSet
from sapo.flowpipe import FlowPipePlotter
from models.sir import SIR

from sapo.benchmark import Benchmark

def test_SIR():

    model = SIR()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(300)

    FlowPipePlotter(mod_flow).plot2DProj(0,1,2)

    print("Average transform time: {0}".format(Benchmark.avg_transf_time()))
