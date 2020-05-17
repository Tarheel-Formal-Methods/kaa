from sapo.reach import ReachSet
from sapo.flowpipe import FlowPipePlotter
from models.phos import Phosphorelay

def test_Phos():

    model = Phosphorelay()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(50)

    FlowPipePlotter(mod_flow).plot2DProj(0,1,2)
