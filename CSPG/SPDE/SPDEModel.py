import numpy as np
from progressbar import Bar, ETA, Percentage, ProgressBar

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "bouchot@mathc.rwth-aachen.de"
__status__ = "Development"
__lastmodified__ = "2015/09/21"

class SPDEModel:
    def samples(self, Z, ratio=None):

        if ratio is not None:
            self.refine_mesh(ratio)

        m, _ = Z.shape

        # Show a progressbar
        widgets = [Percentage(), ' ', Bar(), ' ', ETA()]
        pbar    = ProgressBar(widgets=widgets)

        # Compute all functional evalutations with given precision
        y = np.zeros(m)
        for k in pbar(range(m)):
            y[k] = self.sample(Z[k])

        return y
