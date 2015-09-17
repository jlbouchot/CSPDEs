from dolfin import *

class AbsSquared:
    # Returns form calculating average of the solution for given model
    def __call__(self, model, u, measure):
        return abs(u)**2 * measure
