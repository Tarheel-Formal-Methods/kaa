from kaa.reach import ReachSet
from kaa.flowpipe import FlowPipePlotter

from models.rossler import Rossler

def test_phase_rossler():

    basic_mod = Rossler()
    basic_reach = ReachSet(basic_mod)

    flowpipe = basic_reach.computeReachSet(100)
    FlowPipePlotter(flowpipe).plot2DPhase(0,1)
