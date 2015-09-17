import numpy as np

# Name derived from conditions in Section 6 of "Sparse Legendre Expansions via l1_minimization"
class LD_bounded_operator:
    # TODO: generalize to arbitrary intervals ?
    @staticmethod
    def apply_precondition_measure(v):
        return np.cos(np.pi * v)

