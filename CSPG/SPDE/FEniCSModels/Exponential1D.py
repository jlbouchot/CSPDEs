from dolfin import *

class Exponential1D:
    # TODO: pass min_interval, max_interval on call
    def __init__(self, focus, sigma, min_interval, max_interval):
        self.focus = focus
        self.sigma = sigma

        self.min_interval = min_interval
        self.max_interval = max_interval


    def __call__(self, model, u, measure):
        return exp(-abs(SpatialCoordinate(model.mesh)[0]-self.focus)/self.sigma) / (self.sigma * (2. - exp(-(self.max_interval - self.focus)/self.sigma) - exp((self.min_interval - self.focus)/self.sigma))) * u * measure
