import time
from functools import reduce
from operator import add

"""
Timer object containing duration data.
"""
class TimerData:

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
        
"""
Static class containing all timing utility functions and statistics generating routines.
"""
class Timer:

    timer_stack = []
    time_table = {}

    @staticmethod
    def start(label):
        
        if not Timer.timer_stack or Timer.timer_stack[-1].label != label:
            start_timer = TimerData(label)
            start_timer.start()
            Timer.timer_stack.append(start_timer)
        else:
            raise RuntimeError("Timer of same category started previously for Timer: {}.".format(label))

    @staticmethod
    def stop(label):

        if Timer.timer_stack[-1].label == label:
            end_timer = Timer.timer_stack.pop()
            end_timer.end()

            if end_timer.label not in Timer.time_table:
                Timer.time_table[end_timer.label] = []

            Timer.time_table[end_timer.label].append(end_timer)
            return end_timer.duration
        else:
            raise RuntimeError("Previous timer has not been stopped yet or timer has not been instantiated for Timer: {}.".format(label))

    @staticmethod
    def generate_stats():

        for label, times in Timer.time_table.items():
            print("Average {} Duration: {} sec".format(label, Timer.avg_time(times)))

    def avg_time(times):
        return reduce(add, [t.duration for t in times]) / len(times)


    def total_time(times):
        return reduce(add,[t.duration for t in times])
