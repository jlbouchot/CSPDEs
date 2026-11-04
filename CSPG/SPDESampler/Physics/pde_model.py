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

from dolfin import Constant, TestFunction, TrialFunction, dx, grad, inner

class Physics:
    '''UFL forms of the problem. Declares its structure; the sampler is chosen from it.'''
    num_params = 0
    affine = False            # coefficient affine in z (verified numerically by AffineLinearSampler if needed)
    nonlinear = False
    symmetric = True          # (Jacobian) matrix symmetric positive definite on the free dofs => Cholesky allowed

    def coefficient(self, x, z):
        raise NotImplementedError

    def bilinear(self, disc, c):            # linear problems: a(w, v)
        raise NotImplementedError

    def linear(self, disc):                 # right-hand side L(v)
        raise NotImplementedError

    def residual(self, disc, u, c):         # non-linear problems: F(u; v)
        raise NotImplementedError


class Diffusion(Physics):
    '''-div(a(x,z) grad u) = f, the stereotypical linear parametric problem. Its Jacobian is symmetric positive definite, so Cholesky solver can be used.'''
    affine = True

    def __init__(self, coefficient, f=10.0):
        self.coef, self.f = coefficient, f
        self.num_params = coefficient.num_params

    def coefficient(self, x, z):
        return self.coef(x, z)

    def bilinear(self, disc, c):
        return c * inner(grad(TrialFunction(disc.V)), grad(TestFunction(disc.V))) * dx

    def linear(self, disc):
        return Constant(self.f) * TestFunction(disc.V) * dx

    def residual(self, disc, u, c):         # TODO: lets the Newton sampler solve linear problems too (1 step)
        v = TestFunction(disc.V)
        return c * inner(grad(u), grad(v)) * dx - Constant(self.f) * v * dx


class NonlinearDiffusion(Diffusion):
    '''
    -div((a(x,z) + kappa u^2) grad u) = f: example of a non-linear parametric problem.
    Its Jacobian has the term 2 kappa u du grad(u).grad(v): NOT symmetric, so no Cholesky solver.
    '''
    affine, nonlinear, symmetric = False, True, False

    def __init__(self, coefficient, f=10.0, kappa=1.0):
        super().__init__(coefficient, f)
        self.kappa = kappa

    def residual(self, disc, u, c):
        v = TestFunction(disc.V)
        return (c + Constant(self.kappa) * u ** 2) * inner(grad(u), grad(v)) * dx - Constant(self.f) * v * dx
