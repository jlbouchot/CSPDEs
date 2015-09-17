import numpy as np
from progressbar import Bar, ETA, Percentage, ProgressBar

class SPDEModel:
    def samples(self, Z, ratio=None):

        if ratio is not None:
            self.refine_mesh()

        m, _ = Z.shape

        # Show a progressbar
        widgets = [Percentage(), ' ', Bar(), ' ', ETA()]
        pbar    = ProgressBar(widgets=widgets)

        # Compute all functional evalutations with given precision
        y = np.zeros(m)
        for k in pbar(range(m)):
            y[k] = self.sample(Z[k])

        return y
