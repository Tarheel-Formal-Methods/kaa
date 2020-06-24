from kaa.reach import ReachSet
from kaa.flowpipe import FlowPipePlotter
from models.vanderpol import VanDerPol

def test_VDP():

    model = VanDerPol()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(200)

    FlowPipePlotter(mod_flow).plot2DPhase(0,1)
