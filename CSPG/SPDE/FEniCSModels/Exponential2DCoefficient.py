from dolfin import *

class Exponential2DCoefficient:
    # TODO: pass min_interval, max_interval on call
    def __init__(self, focus, sigma):
        self.focus = focus
        self.sigma = sigma

    def __call__(self, x, z):
        return 10*exp(-(pow(x[0] - self.focus[0], 2) + pow(x[1] - self.focus[1], 2)) / self.sigma)
