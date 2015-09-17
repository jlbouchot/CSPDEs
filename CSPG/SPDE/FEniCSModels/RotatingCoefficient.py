import numpy as np

class RotatingCoefficient:
    def __init__(self, c):
        self.c          = c
        self.num_params = len(c)

    def __call__(self, x, z):
        return self.c * (z+1)/2.
