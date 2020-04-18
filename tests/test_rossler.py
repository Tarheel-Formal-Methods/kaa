from sapo.reach import ReachSet
from models.rossler import Rossler

def test_Rossler():

    model = Rossler()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(100)

    
