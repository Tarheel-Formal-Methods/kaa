from kaa.reach import ReachSet
from kaa.flowpipe import FlowPipePlotter

from models.basic.basic import Basic
import kaa.benchmark as Benchmark

def test_plot_basic():

    basic_mod = Basic()
    basic_reach = ReachSet(basic_mod)

    flowpipe = basic_reach.computeReachSet(10)
    FlowPipePlotter(flowpipe).plot2DProj(0)

    Benchmark.generate_stats()
