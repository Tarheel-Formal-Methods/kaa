from sapo.reach import ReachSet
from sapo.flowpipe import FlowPipePlotter

from models.basic import Basic

def test_plot_basic():

    basic_mod = Basic()
    basic_reach = ReachSet(basic_mod)

    flowpipe = basic_reach.computeReachSet(100)
    FlowPipePlotter(flowpipe).plot2DPhase(0,1)
