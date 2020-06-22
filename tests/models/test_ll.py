from sapo.reach import ReachSet
from sapo.flowpipe import FlowPipePlotter
from models.LL import LL

import sapo.benchmark as Benchmark
from sapo.benchmark import Label

def test_LL():

    model = LL()
    mod_reach = ReachSet(model)

    timer = Benchmark.assign_timer(Label.TOTAL)
    timer.start()
    mod_flow = mod_reach.computeReachSet(100)
    timer.end()

    FlowPipePlotter(mod_flow).plot2DProj(3)

    Benchmark.generate_stats()
