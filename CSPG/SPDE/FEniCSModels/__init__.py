__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "bouchot@mathc.rwth-aachen.de"
__status__ = "Development"
__lastmodified__ = "2016/10/07"

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
#     'UncertainBoundariesDiffusionFEMML',
    'PiecewiseConstantDiffusionFEMModelML',
    'SmallPiecewiseConstantDiffusionFEMModelML',
    'Dim5PiecewiseConstantDiffusionFEMModelML'
    'PiecewiseConstantDiffusionFEMModelML2',
    'PiecewiseConstantDiffusionFEMModelML2_small',
    'PiecewiseConstantDiffusionFEMModelML2_large',
    'PiecewiseConstantDiffusionFEMModelML2_tiny',
    'FinFEMModel',
    'HelmholtzFEMModel',

    'ConstantCoefficient',
    'RotatingCoefficient',
#     'BoundaryPiecewiseConstant',
    'PiecewiseConstantCoefficient',
    'Exponential2DCoefficient',
    'PartitionUnityConstantCoefficient',
    'LinearCoefficient',
    'TrigCoefficient', 
    'CosineCoef1D',
    'CosineCoef2D',
    'CosineCoef3D', 
    'WeightedCosine1D',
    'WeightedCosine2D',
    'WeightedCosine3D'
]

# Functionals
from .Average import *
from .Absolute import *
from .AbsSquared import *
from .Integration import *
from .Exponential1D import *
from .ExponentialUnnormalized import *

# Models
from .ConvectionDiffusionFEMModel import *
from .DiffusionFEMModel import *
from .FinFEMModel import *
from .HelmholtzFEMModel import *

# ML Models
from .DiffusionFEMModelML import *
# from UncertainBoundariesDiffusionFEMML import *
from .PiecewiseConstantDiffusionFEMModelML import *
from .SmallPiecewiseConstantDiffusionFEMModelML import *
from .Dim5PiecewiseConstantDiffusionFEMModelML import *
from .PiecewiseConstantDiffusionFEMModelML2D import *
from .PiecewiseConstantDiffusionFEMModelML2D_small import *
from .PiecewiseConstantDiffusionFEMModelML2D_large import *
from .PiecewiseConstantDiffusionFEMModelML2D_tiny import *

# Coefficients
from .ConstantCoefficient import *
from .RotatingCoefficient import *
from .LinearCoefficient import *
from .PiecewiseConstantCoefficient import *
# from BoundaryPiecewiseConstant import *
from .Exponential2DCoefficient import *
from .TrigCoefficient import *
from .CosineCoef1D import *
from .CosineCoef2D import *
from .CosineCoef3D import *
from .PartitionUnityConstantCoefficient import *
from .WeightedCosine1D import *
from .WeightedCosine2D import *
from .WeightedCosine3D import *
