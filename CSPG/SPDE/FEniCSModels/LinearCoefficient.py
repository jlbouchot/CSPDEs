class LinearCoefficient:
    num_params = 1

    def __init__(self, a, b):
        self.a = a
        self.b = b

    def __call__(self, x, z):
        return self.a + (z[0]+1)/2 * (self.b - self.a)
