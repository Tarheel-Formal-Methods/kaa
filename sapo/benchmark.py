import time
from enum import Enum
from functools import reduce
from operator import add

class Label(Enum):
    TOTAL = 1
    TRANSF = 2
    PLOT_PROJ = 3
    PLOT_PHASE = 4

class Timer:

    def __init__(self, label):

        self.label = label
        self.start_time = 0
        self.end_time = 0
        self.duration = 0

    def start(self):
        self.start_time = time.time()

    def end(self):
        self.end_time = time.time()
        self.duration = self.end_time - self.start_time


class Benchmark:

    @staticmethod
    def assign_timer(label):
        new_timer = Timer(label)
        BenchmarkUtils.timer_bank.append(new_timer)

        return new_timer

    @staticmethod
    def avg_transf_time():

        durations = [ timer.duration for timer in BenchmarkUtils.timer_bank if timer.label == Label.TRANSF ]
        avg_duration = (reduce(add, durations)) / len(durations)

        return avg_durations

class BenchmarkUtils:

    timer_bank = []
