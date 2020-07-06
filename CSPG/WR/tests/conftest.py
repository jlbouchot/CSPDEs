import pytest
import numpy as np
import WR
from scipy.sparse import rand

parameter_set = {"epsilon": 1e-6,
                 "density": 0.15,
                 "n": 1000,
                 "maxiter": 50}

# randstate = 42
# m = int(np.log(n) * density * n)
nb_test_runs = 5

def sparse_ground_truth(size, density, randstate):
    matrix = rand(size, 1, density=density, format="csr", random_state=randstate).todense()
    return np.squeeze(np.asarray(matrix)) #squeeze needed because rand returns the deprecated matrix and we want to reduce dim

def create_Operator(m,n,randstate):
    np.random.seed(randstate)
    matrix = np.random.normal(size=(m,n))
    matrix = matrix/(np.sqrt(m))
    operator = WR.Operators.operator_from_matrix(None, matrix)
    return operator

def create_obs_y(op,x):
    return op.apply(x)

@pytest.fixture(params=np.asarray(range(nb_test_runs))+1)
def para_set(request):
    parameter_set["randstate"] = request.param
    parameter_set["m"] = int(np.log(parameter_set["n"]) * parameter_set["density"] * parameter_set["n"])
    parameter_set["x_truth"] = sparse_ground_truth(parameter_set["n"],
                                                   parameter_set["density"],
                                                   parameter_set["randstate"])
    parameter_set["operator"] = create_Operator(parameter_set["m"],
                                                parameter_set["n"],
                                                parameter_set["randstate"])
    parameter_set["y"] = create_obs_y(parameter_set["operator"], parameter_set["x_truth"])
    return parameter_set

@pytest.fixture(params=[WR.whtp, WR.exact_wbp_cvx, WR.qc_wbp_cvx,
                        WR.Qc_wbp_precond_primaldual, WR.wcosamp,
                        WR. wiht, WR.womp, WR.primaldual])
def algo(request):
    return request.param
