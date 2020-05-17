from sapo.reach import ReachSet
from sapo.flowpipe import FlowPipePlotter
from models.vanderpol import VanDerPol

def test_VDP():

    model = VanDerPol()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(200)

    FlowPipePlotter(mod_flow).plot2DPhase(0,1)


if __name__ == '__main__':
    test_VDP()
