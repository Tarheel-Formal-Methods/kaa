from kaa.reach import ReachSet
from kaa.flowpipe import FlowPipePlotter

from models.basic.basic import Basic

def test_phase_basic():

    basic_mod = Basic()
    basic_reach = ReachSet(basic_mod)

    flowpipe = basic_reach.computeReachSet(100)
    FlowPipePlotter(flowpipe).plot2DPhase(0,1)
