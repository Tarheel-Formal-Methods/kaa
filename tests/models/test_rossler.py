from kaa.reach import ReachSet
from kaa.flowpipe import FlowPipePlotter
from models.rossler import Rossler

import kaa.benchmark as Benchmark
from kaa.benchmark import Label

def test_Rossler():

    model = Rossler()
    mod_reach = ReachSet(model)

    timer = Benchmark.assign_timer(Label.TOTAL)

    timer.start()
    mod_flow = mod_reach.computeReachSet(300)
    timer.end()

    FlowPipePlotter(mod_flow).plot2DProj(0)
    Benchmark.generate_stats()
