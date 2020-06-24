from kaa.reach import ReachSet
from kaa.flowpipe import FlowPipePlotter
from models.lotkavolterra import LotkaVolterra

import kaa.benchmark as Benchmark
from kaa.benchmark import Label

def test_LV():

    model = LotkaVolterra()
    mod_reach = ReachSet(model)

    timer = Benchmark.assign_timer(Label.TOTAL)
    timer.start()

    mod_flow = mod_reach.computeReachSet(300)

    timer.end()

    FlowPipePlotter(mod_flow).plot2DProj(0)
    FlowPipePlotter(mod_flow).plot2DProj(1)
    FlowPipePlotter(mod_flow).plot2DProj(2)
    Benchmark.generate_stats()
