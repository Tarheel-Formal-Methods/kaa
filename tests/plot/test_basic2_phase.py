from kaa.reach import ReachSet
from kaa.flowpipe import FlowPipePlotter

from models.basic.basic2 import Basic2

def test_phase_basic2():

    basic_mod = Basic2()
    basic_reach = ReachSet(basic_mod)

    flowpipe = basic_reach.computeReachSet(100)
    FlowPipePlotter(flowpipe).plot2DPhase(0,1)
