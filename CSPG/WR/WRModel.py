import numpy as np

class WRModel:
    def __init__(self, method, operator, weights, m_from_s_N, check):
        self.method         = method
        self.operator       = operator
        self.weights        = weights
        self.get_m_from_s_N = m_from_s_N
        self.check          = check

    def estimate_ML_samples(self, cspde_result, Z_cross):
        nbLevel = len(cspde_result)
        y_recon = np.zeros(len(Z_cross))

        for oneLevel in xrange(0,nbLevel):
            # print cspde_result[oneLevel].result.x.shape
            # print cspde_result[oneLevel].result.x[cspde_result[oneLevel].result.x != 0]
            # print cspde_result[oneLevel].J_s.shape
            F_recon = self.operator.create(np.array(cspde_result[oneLevel].J_s)[cspde_result[oneLevel].result.x != 0], Z_cross, normalization=np.sqrt(cspde_result[oneLevel].m))
            y_recon = y_recon+F_recon.apply(cspde_result[oneLevel].result.x[cspde_result[oneLevel].result.x != 0])

        return y_recon


def check_cs(N, m):
    assert N >= m, "It is N < m ({0} < {1}). Aborting ...".format(N, m)

def cs_theoretic_m(s, N):
    # TODO: Consant factor is known to be 2?
    return int(np.ceil(2 * np.log(s)**3 * s * np.log(N)))

def cs_pragmatic_m(s, N):
    # Calculate number of samples using the "pragmatic" bound. That is:
    #  Factor np.log(s)**3 probably ignoreable in practice
    #  Constant factor "C" is about 2
    return int(np.ceil(2 * s * np.log(N)))
	
def cs_theoretic_m_new(s, N):
    # TODO: Consant factor is known to be 2?
    return int(np.ceil(2 * np.log(s)**2 * s * np.log(N)))