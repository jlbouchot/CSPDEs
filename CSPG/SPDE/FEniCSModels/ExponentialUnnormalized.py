from dolfin import *

import numpy as np

class ExponentialUnnormalized:
    # TODO: pass min_interval, max_interval on call
    def __init__(self, focus, sigma):
        self.focus = focus
        self.sigma = sigma

    def __call__(self, model, u):
        x = SpatialCoordinate(model.mesh)

        return exp(-abs(sum([x[k] - self.focus[k] for k in range(model.mesh.topology().dim())]))/self.sigma) * u * dx(model.mesh)
