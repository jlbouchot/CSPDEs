from test_diff_avg_v_ML_2D import Main as DIFF_2D

# Global parameters: 
uniform_weights = 1.08
sampling_name = "p"
nb_iter = 50
epsilon = 1e-4
nb_tests = 10

alpha = 2.0
abar = 4.3

algos = ['whtp']
dimensions = [10]
L = 1
dat_cs = [5,7,10,14,20,28,40,56,80]

n0s = [28] # A reasonably coarse / fine grid

for recovery in algos:
    print("\n \t***** Let's start with {} as a recovery algorithm *****".format(recovery))
    for oneD in dimensions:
        print("\n\t\t***** We'll be working in {} dimensions".format(oneD))
        for n0 in n0s:
            for dat_constant in dat_cs:
                # Generate meaningful output file name
                outputFile = '_'.join(['diffusion', 'cosines', '2D', 'd', str(oneD), 'n0', str(n0), 'c', str(dat_constant), 'v', str(uniform_weights), 'L', str(L), recovery])
                # 
                DIFF_2D(outputFile, oneD, tuple([n0, n0]), L, recovery, uniform_weights, L, sampling_name, nb_iter, epsilon, nb_tests, alpha, abar, dat_constant)
