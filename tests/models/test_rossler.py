from kaa.reach import ReachSet
from kaa.flowpipe import FlowPipePlotter
from models.rossler import Rossler

from kaa.timer import Timer

def test_Rossler():

    model = Rossler()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(200)

    Timer.generate_stats()
