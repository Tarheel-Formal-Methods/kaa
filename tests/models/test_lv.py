from sapo.reach import ReachSet
from sapo.flowpipe import FlowPipePlotter
from models.lotkavolterra import LotkaVolterra

import sapo.benchmark as Benchmark

def test_LV():

    model = LotkaVolterra()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(100)

    FlowPipePlotter(mod_flow).plot2DPhase(0,1)
    Benchmark.generate_stats()
