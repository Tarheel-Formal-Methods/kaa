from sapo.reach import ReachSet
from sapo.flowpipe import FlowPipePlotter
from models.sir import SIR

import sapo.benchmark as Benchmark

def test_SIR():

    model = SIR()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(300)

    FlowPipePlotter(mod_flow).plot2DProj(0,1,2)

    Benchmark.generate_stats()
