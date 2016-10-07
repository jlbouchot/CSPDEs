Here is a list of things still to be done (with various level of priorities and difficulties)
=============================================================================================

Code organization: 
------------------
* Separate tests from actual developments (include all the test_XXX files)
* Have some fairly simple examples

Code improvements:
------------------
* Automatize the piecewise linear diffusion problem to have a dynamic number of splits
* Allow for uncertainty in where these splits are actually
* Implement some parallelization of the calculations of the samples
* Implement a tensor-based calculation of the matrix (i.e. avoid the actual computations / storage of the matrix)

New ideas: 
----------
* Extend to time-dependent equations
* Embed some information about the regularity (expected) of the solution
* What about non-independent draws of samples from one level to another? (this would save a little bit of computation time)

General things:
---------------
* Write a simple doc and hands-on / guide for further uses
