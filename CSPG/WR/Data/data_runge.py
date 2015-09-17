class data_runge:
    name = 'Runge Function'

    def __init__(self, t):
        self.t = t
        self.y = 1.0 / (1 + 25 * t ** 2)
