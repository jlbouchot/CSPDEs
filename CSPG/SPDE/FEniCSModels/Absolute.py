from dolfin import *

class Absolute:
    # Returns form calculating average of the solution for given model
    def __call__(self, model, u, measure):
        return abs(u) * measure
