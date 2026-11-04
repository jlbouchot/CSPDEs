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

from abc import ABC, abstractmethod


class AbstractQOI(ABC):
    """Abstract base class for quantities of interest."""

    @abstractmethod
    def form(self, disc, u):
        """Return the variational form of this QoI for ``u`` on ``disc``, a discretisation of the PDE."""
        raise NotImplementedError
