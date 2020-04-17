from sapo.reach import ReachSet
from sapo.sir import SIR

def main():

    model = SIR()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(300)


if __name__ == "__main__":
    main()
