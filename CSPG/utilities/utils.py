# A utility package for all sorts of helper functions for MLCSPG
import WR

def get_sampling_type(sampling_name):
    switcher = {
        "pragmatic": WR.cs_pragmatic_m,
        "theoretic": WR.cs_theoretic_m,
        "new": WR.cs_theoretic_m_new,
		"p": WR.cs_pragmatic_m,
        "t": WR.cs_theoretic_m,
    }
    return switcher.get(sampling_name, WR.cs_pragmatic_m)

def is_tensor_based_possible(method):
    """ 
    Returns True if the method is compatible with tensor-based operators, False otherwise. 
    It is false for anything based on convex optimization, as these methods are not compatible with tensor-based operators.
    """
    return method not in ["bp", "bpdn"]