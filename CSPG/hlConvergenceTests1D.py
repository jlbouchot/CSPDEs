from test_diff_avg_v_ML import Main as DIFF_1D
from test_diff_avg_v_ML_2D import Main as DIFF_2D
from test_diff_avg_v_ML_3D import Main as DIFF_3D


# Global parameters: 
L = 3
minL = 3 # Test only one case
dat_constant = 8
# recovery = "bp"
uniform_weights = 1.08
sampling_name = "p"
nb_iter = 50
epsilon = 1e-4
nb_tests = 10

alpha = 2.0
abar = 4.3

# what to test: 
n0s = [5,7,10,14,20,28,40,56,70]
dimensions = [10,15,20] # Will probably split into  different files to avoid taking too much time -- this shouldn't be a problem for 2 and 1 D cases. 
algos = ['bp', 'womp', 'whtp']

for recovery in algos:
    print("\n \t***** Let's start with {} as a recovery algorithm *****".format(recovery))
    for oneD in dimensions:
        print("\n\t\t***** We'll be working in {} dimensions".format(oneD))
        for n0 in n0s:
            # Generate meaningful output file name
            outputFile = '_'.join(['diffusion', 'cosines', '1D', 'd', str(oneD), 'n0', str(n0), 'c', str(dat_constant), 'v', str(uniform_weights), 'L', str(L), recovery])
            # 
            DIFF_1D(outputFile, oneD, tuple([n0]), L, recovery, uniform_weights, L, sampling_name, nb_iter, epsilon, nb_tests, alpha, abar, dat_constant)


