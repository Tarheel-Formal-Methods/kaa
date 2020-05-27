from sapo.reach import ReachSet
from sapo.flowpipe import FlowPipePlotter
from models.phos import Phosphorelay

import sapo.benchmark as Benchmark
from sapo.benchmark import Label

def test_Phos():

    model = Phosphorelay()
    mod_reach = ReachSet(model)
    timer = Benchmark.assign_timer(Label.TOTAL)
    timer.start()

    mod_flow = mod_reach.computeReachSet(300)

    timer.end()

    FlowPipePlotter(mod_flow).plot2DProj(0,1,2)
    Benchmark.generate_stats()
