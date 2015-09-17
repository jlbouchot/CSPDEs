import numpy as np

class Constant:
    def __init__(self, c):
        self.name = "all-{0}".format(c)
        self.w    = c * np.ones(100000)
