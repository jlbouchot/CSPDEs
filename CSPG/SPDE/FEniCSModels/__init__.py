__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "bouchot@mathc.rwth-aachen.de"
__status__ = "Development"
__lastmodified__ = "2015/09/22"

__all__ = [
    'Average',
    'Absolute',
    'AbsSquared',
    'Integration',
    'Exponential1D',
    'ExponentialUnnormalized',

    'ConvectionDiffusionFEMModel',
    'DiffusionFEMModel',
    'DiffusionFEMModelML',
    'PiecewiseConstantDiffusionFEMModelML'
    'FinFEMModel',
    'HelmholtzFEMModel',

    'ConstantCoefficient',
    'RotatingCoefficient',
    'PiecewiseConstantCoefficient',
    'Exponential2DCoefficient',
    'PartitionUnityConstantCoefficient',
    'LinearCoefficient',
    'TrigCoefficient'
]

# Functionals
from Average import *
from Absolute import *
from AbsSquared import *
from Integration import *
from Exponential1D import *
from ExponentialUnnormalized import *

# Models
from ConvectionDiffusionFEMModel import *
from DiffusionFEMModel import *
from FinFEMModel import *
from HelmholtzFEMModel import *

# ML Models
from DiffusionFEMModelML import *
from PiecewiseConstantDiffusionFEMModelML import *

# Coefficients
from ConstantCoefficient import *
from RotatingCoefficient import *
from LinearCoefficient import *
from PiecewiseConstantCoefficient import *
from Exponential2DCoefficient import *
from TrigCoefficient import *
from PartitionUnityConstantCoefficient import *
