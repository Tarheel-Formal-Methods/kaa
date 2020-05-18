from sapo.reach import ReachSet
from sapo.flowpipe import FlowPipePlotter

from models.sir import SIR

def test_phase_sir():

    sir_mod = SIR()
    sir_reach = ReachSet(sir_mod)

    flowpipe = sir_reach.computeReachSet(100)
    FlowPipePlotter(flowpipe).plot2DPhase(0,1)
