import numpy as np

class data_sparse_power_distributed:
    def __init__(self, Operator, s, alpha):
        self.name = 'Simple Power Law Data from ' + Operator.name

        x      = np.zeros(Operator.n, 1)
        x[0:s] = np.arange(1,s+1)**alpha
        self.x = x

        self.y = Operator.apply(x)
