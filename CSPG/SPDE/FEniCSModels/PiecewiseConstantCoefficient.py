from ufl import *

class PiecewiseConstantCoefficient:
    def __init__(self, c):
        self.num_params = c

    def __call__(self, x, z):
        def cond(k):
            if k >= len(self.c):
                return self.c[-1][0]

            return conditional(ge(self.c[k][0], x[0]), self.c[k][1], cond(k+1))

        return cond(0)
