class ConstantCoefficient:
    num_params = 0

    def __init__(self, c):
        self.c = c

    def __call__(self, x, z):
        return self.c
