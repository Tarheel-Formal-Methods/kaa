from sapo.reach import ReachSet
from sapo.flowpipe import FlowPipePlotter
from models.vanderpol import VanDerPol

def test_VDP():

    model = VanDerPol()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(100)

    FlowPipePlotter(mod_flow).plot2DProj(0)
