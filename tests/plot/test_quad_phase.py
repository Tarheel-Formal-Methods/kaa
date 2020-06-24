from kaa.reach import ReachSet
from kaa.flowpipe import FlowPipePlotter

from models.quadcopter import Quadcopter

def test_phase_quad():

    quad_mod = Quadcopter()
    quad_reach = ReachSet(quad_mod)

    flowpipe = quad_reach.computeReachSet(50)
    FlowPipePlotter(flowpipe).plot2DPhase(1,3)
