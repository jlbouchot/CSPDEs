import numpy as np
import os.path

from .operator_from_matrix_Alt import *
from .LD_bounded_operator  import LD_bounded_operator

__author__ = ["Falk Pulsmeyer", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2026, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Falk Pulsmeyer", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__lastmodified__ = "2026/07/21"

class Cheb_Alt(LD_bounded_operator):
    # L_inf norm of the basis functions associated to this operator
    theta = np.sqrt(2)
    name  = 'Cheb_Alt'

    @staticmethod
    def create(J, Z, normalization=None):
        def base(x, k):
            return np.cos(k * np.arccos(x)) * np.sqrt(2)**(k>0)
        return operator_from_matrix_Alt(Cheb_Alt, matrix_from_tensor_indices(J, Z, base, normalization), univ_tensor_from_tensor_indices(J, Z, base, normalization),J)

    @staticmethod
    def load(data_file, t_fname):
        dirname = os.path.dirname(data_file)
        bname = os.path.basename(data_file)
        A = np.load(os.path.join(dirname, "A_"+bname))
        # Univariate blocks were saved as "{idx}_univariate_<bname>.npy". Collect and load them.
        univariate_files = [f for f in os.listdir(dirname) if f.endswith("_univariate_"+bname)]
        # Sort by the leading index to preserve original order (assumes integer prefixes)
        univariate_files.sort(key=lambda s: int(s.split("_univariate_")[0]))
        univariate = []
        for fname in univariate_files:
            univariate.append(np.load(os.path.join(dirname, fname)))
        # univariate = np.load(os.path.join(dirname, "univariate_"+bname))
        J = np.load(os.path.join(dirname, "J_"+bname))

        return operator_from_matrix_Alt(Cheb_Alt, A, univariate, J), np.load(t_fname)
