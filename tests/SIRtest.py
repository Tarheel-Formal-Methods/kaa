from sapo.reach import ReachSet
from models.sir import SIR

def test_SIR():

    model = SIR()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(300)
