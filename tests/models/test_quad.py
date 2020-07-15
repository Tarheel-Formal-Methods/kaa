from kaa.reach import ReachSet
from kaa.flowpipe import FlowPipePlotter
from models.quadcopter import Quadcopter

from kaa.timer import Timer

def test_Quad():

    model = Quadcopter()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(300)

    FlowPipePlotter(mod_flow).plot2DProj(2)
    FlowPipePlotter(mod_flow).plot2DProj(5)
    FlowPipePlotter(mod_flow).plot2DProj(13)
    Timer.generate_stats()
