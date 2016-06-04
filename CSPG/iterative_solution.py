import numpy as np

def compute_true_avg(ys, rhs, variability, mean_field):
    
    nb_param = 13 # We are dealing with 13 pieces
    delta_ts = 1.0/(nb_param) # length of the segments
    print "Delta = ", delta_ts
    if type(ys) is list:
        ys = np.array(ys) 
    if len(ys.shape) == 2:
        dim,nb_ys = ys.shape
    else: 
        dim = ys.shape[0]
        nb_ys = 1

    one_over_a = np.divide(1.0,mean_field+variability*ys)
    harm_means = np.cumsum(one_over_a,axis = 0)
    weighted_means = np.cumsum(np.array([2*np.linspace(1,nb_param,nb_param)-1]*nb_ys).transpose()*one_over_a,axis = 0)

    C_2 = np.divide(1.0,harm_means[-1,:])/delta_ts + rhs*delta_ts / 2.0 * np.divide(weighted_means[-1,:],harm_means[-1,:])

    u_tilda_k = delta_ts*np.array([C_2]*nb_param)*one_over_a - rhs*delta_ts**2/2*np.array([2*np.linspace(1,nb_param,nb_param)-1]*nb_ys).transpose()*one_over_a
    u_tilda_k = np.vstack((np.zeros((1,nb_ys), dtype=float), u_tilda_k ) )
    u_k = np.cumsum(u_tilda_k, axis = 0)
    to_add = delta_ts**2/2*np.array([C_2]*nb_param)*one_over_a - rhs/6.0*delta_ts**3*one_over_a*np.array([3.0*np.linspace(1,nb_param,nb_param)-2.0]*nb_ys).transpose()

    return np.sum(delta_ts*u_k[0:-1,:] + to_add, axis = 0)
