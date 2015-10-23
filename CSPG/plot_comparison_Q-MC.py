from dolfin import *
import ghalton

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import sys
import shelve


def plot_qmc(infile, outfile, v):
    def average_mc(y):
        return y.sum() / len(y)

    print("Loading {0} ...".format(infile))
    results = sorted([r for r in shelve.open(infile).values() if r.wr_model.weights[0] == v], key=lambda r: r.s)

    first_result = results[0]
    d            = first_result.cspde_result.d
    spde_model   = first_result.spde_model
    epsilon      = first_result.epsilon

    # Compute samples
    sequencer = ghalton.Halton(d)
    x_m       = [r.cspde_result.m for r in results]
    y_m_g0    = [r.cspde_result.result.x[0]/np.sqrt(r.cspde_result.m) for r in results]

    y, y_m_mc = [], []
    last_x    = 0
    for x in x_m:
        y = np.hstack((y, spde_model.samples(first_result.wr_model.operator.apply_precondition_measure(2*(np.array(sequencer.get(x - last_x))-.5)), epsilon)))
        y_m_mc.append(average_mc(y))

        last_x = x

    plt.plot(x_m, y_m_mc, x_m, y_m_g0)
    plt.savefig(outfile+".svg", bbox_inches="tight")
    plt.clf()

### Main
if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: {0} test_data output_file v".format(sys.argv[0]))

    plot_qmc(sys.argv[1], sys.argv[2], float(sys.argv[3]))
