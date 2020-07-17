from kaa.reach import ReachSet
from kaa.flowpipe import FlowPipePlotter
from models.LL import LL

from kaa.timer import Timer

def test_LL():

    model = LL()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(150)


    FlowPipePlotter(mod_flow).plot2DProj(3)

    Timer.generate_stats()
