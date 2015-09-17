import numpy as np

#import matplotlib.pyplot as plt
import time


class Test:
    def __init__(self, method, operator, data, weights, *varargin):
        self.method                    = method
        self.operator                  = operator
        self.data                      = data
        self.weights                   = weights
        self.further_method_parameters = varargin

    def test(self):
        tic           = time.time()
        self.solution = self.method(self.operator, self.data.y, self.weights.w[0:self.operator.n], *self.further_method_parameters)
        self.time     = time.time() - tic

    def output_statistics(self):
        print('Reconstructed {0} using {1} algorithm, {2} operator and {3} weights: ||A*x - y||_2 = {4} after {5} iterations and {6} seconds'
              .format(self.data.name, self.solution.methodname, self.operator.model.name, self.weights.name,
              np.linalg.norm(self.operator.apply(self.solution.x) - self.data.y), self.solution.iterations, self.time))

        print('')

    # Works only in the one-dimensional case
#    def output_difference_plot(self,T,xmin,xmax,ymin,ymax,dmin,dmax):
#        TOperator = self.operator.model.create(np.arange(self.operator.n), T, normalization=np.sqrt(self.operator.m))
#        F         = TOperator.apply(self.solution.x)
#        Ydata     = self.data.__class__(T)
#
#        # fig = plt.figure()
#        # fig.subplots_adjust(bottom = 0)
#        # fig.subplots_adjust(top = 2)
#        # fig.subplots_adjust(right = 2)
#        # fig.subplots_adjust(left = 0)
#
#        ax = plt.subplot(211)
#        plt.plot(T,Ydata.y,label="Original function")
#        plt.plot(T,F,label="Reconstructed function using {0} operator, {1} algorithm and {2} weights".format(self.operator.model.name,self.solution.methodname,self.weights.name))
#        plt.plot(self.data.t,self.data.y,'o',label="Sample points")
#        plt.axis([xmin,xmax,ymin,ymax])
#        plt.xlabel('x')
#        plt.ylabel('{0}(x)'.format(self.data.name))
#
#        handles, labels = ax.get_legend_handles_labels()
#        ax.legend(handles, labels,
#                  ncol=1, mode="expand", borderaxespad=0.)#bbox_to_anchor=(0., 1.02, 1., .102), 
#
#        plt.subplot(212)
#        plt.plot(T,Ydata.y - F)
#        plt.axis([xmin,xmax,dmin,dmax])
#        plt.title('Difference')
#        plt.xlabel('x')
#        plt.ylabel('{0}(x) - {1}[Reconstructed](x)'.format(self.data.name,self.data.name))
#
#        plt.gcf().set_size_inches(17, 12)
#        plt.tight_layout()
#
#        plt.savefig('test-{0}-{1}-{2}-{3}.pdf'.format(self.solution.methodname,self.data.name,self.operator.model.name,self.weights.name),
#                    format='pdf')
#
#        plt.clf()
