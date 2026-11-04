__author__ = ["Jean-Luc Bouchot"]
__copyright__ = "Copyright 2017-2026, INRIA, LMU Munich, and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology, INRIA Sophia Antipolis"
__credits__ = ["Jean-Luc Bouchot"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__created__ = "2026/10/29"
__lastmodified__ = "2026/11/04"

from dolfin import *
from .abstract_qoi import AbstractQOI


class AverageQOI(AbstractQOI):
    '''Quantity of interest that computes the average of a function over the domain.'''
    linear = True

    def form(self, disc, u):
        '''
        Return the variational form of this QoI for ``u`` on ``disc``, a discretisation of the PDE.
        A flag indicating that this QoI is computed on a unit domain allows slightly faster computation of the average.
        '''
        if not disc.is_unit:
            domain = u.function_space().mesh()
            volume = assemble(Constant(1.0) * dx(domain=domain))
        else: 
            volume = 1.0
        return (u / volume) * dx
