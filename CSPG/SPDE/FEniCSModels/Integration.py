from dolfin import *

class Integration:
    # Returns form calculating average of the solution for given model
    def __call__(self, model, u, measure):
        return u * measure
