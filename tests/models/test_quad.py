from sapo.reach import ReachSet
from sapo.flowpipe import FlowPipePlotter
from models.quadcopter import Quadcopter

import sapo.benchmark as Benchmark
from sapo.benchmark import Label

def test_Quad():

    model = Quadcopter()
    mod_reach = ReachSet(model)

    timer = Benchmark.assign_timer(Label.TOTAL)
    timer.start()

    mod_flow = mod_reach.computeReachSet(300)

    timer.end()

    FlowPipePlotter(mod_flow).plot2DProj(2,5,13)
    Benchmark.generate_stats()
