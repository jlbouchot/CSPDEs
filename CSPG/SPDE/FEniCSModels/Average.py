from dolfin import *

class Average:
    # Returns form calculating average of the solution for given model
    def __call__(self, model, u, measure):
        # Careful: Only works if all elements of the domain are marked 0
        volume = assemble(Constant(1.0) * measure(model.mesh))

        return 1.0/volume * u * measure
