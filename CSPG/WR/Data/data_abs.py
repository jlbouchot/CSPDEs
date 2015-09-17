import numpy as np

class data_abs:
    name = 'Modulus'

    def __init__(self, t):
        self.t = t
        self.y = np.abs(t)
