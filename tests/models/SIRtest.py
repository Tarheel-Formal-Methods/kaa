from kaa.reach import ReachSet
from kaa.flowpipe import FlowPipePlotter
from models.sir import SIR

import kaa.benchmark as Benchmark
from kaa.benchmark import Label

def test_SIR():

    model = SIR()
    mod_reach = ReachSet(model)

    timer = Benchmark.assign_timer(Label.TOTAL)
    timer.start()
    mod_flow = mod_reach.computeReachSet(300)
    timer.end()

    FlowPipePlotter(mod_flow).plot2DProj(3)

    Benchmark.generate_stats()
