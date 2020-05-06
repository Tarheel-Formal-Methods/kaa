from sapo.reach import ReachSet
from sapo.flowpipe import FlowPipePlotter
from models.quadcopter import Quadcopter

def test_Quad():

    model = Quadcopter()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(200)

    FlowPipePlotter(mod_flow).plot2DProj(2)
