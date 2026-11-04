from dolfin import *
from ufl import And, ge, le
from .abstract_qoi import AbstractQOI

__author__ = ["Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015 - 2026, INRIA, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__created__ = "2026/11/02"
__lastmodified__ = "2026/11/04"

class BoxAverage(AbstractQOI):
    '''
    Average of the solution over the box [x0, x1] x [y0, y1] (default: [0, 1/4]^2). 
    A smooth L^2 functional
    Keep the box edges on mesh lines of the coarsest level (multiples of 1/(10 2^J)) for exact quadrature.
    '''
    linear = True

    def __init__(self, x0 = 0.0, x1 = 0.25, y0 = 0.0, y1 = 0.25):
        self.box = (x0, x1, y0, y1)

    def __call__(self, model, u, measure):
        x0, x1, y0, y1 = self.box
        x = SpatialCoordinate(model.mesh)
        chi = conditional(And(And(ge(x[0], x0), le(x[0], x1)), And(ge(x[1], y0), le(x[1], y1))), 1.0, 0.0)
        return chi * u / ((x1 - x0) * (y1 - y0)) * measure

    def form(self, disc, u):
        x0, x1, y0, y1 = self.box
        x = SpatialCoordinate(disc.mesh)
        chi = conditional(And(And(ge(x[0], x0), le(x[0], x1)), And(ge(x[1], y0), le(x[1], y1))), 1.0, 0.0)
        return chi*u / ((x1 - x0) * (y1 - y0)) * dx