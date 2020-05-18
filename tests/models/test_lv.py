from sapo.reach import ReachSet
from sapo.flowpipe import FlowPipePlotter
from models.lotkavolterra import LotkaVolterra

def test_LV():

    model = LotkaVolterra()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(50)

    FlowPipePlotter(mod_flow).plot2DProj(0,1)
