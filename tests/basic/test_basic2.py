from sapo.reach import ReachSet
from sapo.flowpipe import FlowPipePlotter

from models.basic2 import Basic2

def test_basic2():

    basic_mod = Basic2()
    basic_reach = ReachSet(basic_mod)

    flowpipe = basic_reach.computeReachSet(150)
    FlowPipePlotter(flowpipe).plot2DProj(0)
