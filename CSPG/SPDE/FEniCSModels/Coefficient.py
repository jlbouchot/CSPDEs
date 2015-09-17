class Coefficient:
    num_params = 0
    # @ lru_cache(...)
    def __call__(self, x, z):
        return self.a_bar(x) + self.phi(x, z)
