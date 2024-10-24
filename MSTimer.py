from MSLogging import Logger
import time
class Timer:
    def __init__(self, timer_name: ""):
        self.timer_name = timer_name
        self.timers = {}
        self.prev = time.time()

    def init(self):
        self.prev = time.time()
        self.timers = {}

    def reset(self):
        self.prev = time.time()

    def add(self, *names):
        for name in names:
            if name in self.timers:
                Logger.warning("Timer: timer already exists: " + name)
                continue
            self.timers[name] = 0.0

    def elapsed_and_reset(self, name):
        if name not in self.timers:
            self.timers[name] = 0.0

        self.timers[name] += time.time() - self.prev
        self.prev = time.time()

    def print(self):
        for name, elapsed in self.timers.items():
            Logger.debug(f"[{self.timer_name}] " + name + ": " + "%.2f" % elapsed + "s")
