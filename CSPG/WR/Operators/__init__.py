__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2019, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Falk Pulsmeyer", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__lastmodified__ = "2026/10/14"

__all__ = [
    'Chebyshev',
    'Cos',
    'Legendre',
    'Cheb_Alt',
    'Cheb_AltPP',
    'Cheb_Sparse'
]

from .Chebyshev import *
from .Cheb_Alt import *
from .Cheb_AltPP import *
from .Cheb_Sparse import *
from .Cos import *
from .Legendre import *
