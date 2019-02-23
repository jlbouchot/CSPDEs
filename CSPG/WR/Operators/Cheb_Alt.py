import numpy as np

from .operator_from_matrix_Alt import *
from .LD_bounded_operator  import LD_bounded_operator

__author__ = ["Falk Pulsmeyer", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2019, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Falk Pulsmeyer", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__lastmodified__ = "2019/02/22"

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
    def load(data_file):
        previous_data = np.load(data_file)

        return operator_from_matrix_Alt(Cheb_Alt, previous_data[0], previous_data[1], previous_data[2])
