import numpy as np

class data_sign:
    name = 'Sign'

    def __init__(self, t):
        self.t = t
        self.y = np.sign(t)
