from kaa.reach import ReachSet
from kaa.flowpipe import FlowPipePlotter
from models.ebola import Ebola

def test_Ebola():
    model = Ebola()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(100)

    FlowPipePlotter(mod_flow).plot2DProj(0,1)
