from kaa.reach import ReachSet
from kaa.flowpipe import FlowPipePlotter
from kaa.timer import Timer

from models.basic.basic2 import Basic2

def test_basic2():

    basic_mod = Basic2()
    basic_reach = ReachSet(basic_mod)

    flowpipe = basic_reach.computeReachSet(300)
    FlowPipePlotter(flowpipe).plot2DProj(0)

    Timer.generate_stats()
