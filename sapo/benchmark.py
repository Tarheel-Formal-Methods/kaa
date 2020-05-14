import time
from enum import Enum, auto
from functools import reduce
from operator import add

class Label(Enum):
    TOTAL = auto()
    TRANSF = auto()
    PLOT_PROJ = auto()
    PLOT_PHASE = auto()
    BERN = auto()

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


timer_bank = []

def assign_timer(label):
    new_timer = Timer(label)
    timer_bank.append(new_timer)

    return new_timer

def generate_stats():
    print('Average Bundle Transformation Duration: {0}'.format(avg_time(Label.TRANSF)))
    print('Average Bernstein Computation Duration: {0}'.format(avg_time(Label.BERN)))
    print('Average Phase Computation Duration: {0}'.format(avg_time(Label.PLOT_PHASE)))

def avg_time(label):

    durations = [ timer.duration for timer in timer_bank if timer.label == label ]
    avg_duration = (reduce(add, durations)) / len(durations)
    return avg_duration
