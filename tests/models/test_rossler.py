from sapo.reach import ReachSet
from sapo.flowpipe import FlowPipePlotter
from models.rossler import Rossler

def test_Rossler():

    model = Rossler()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(100)

    FlowPipePlotter(mod_flow).plot2DProj(2)
