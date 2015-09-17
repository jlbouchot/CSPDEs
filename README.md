CSPDEs
======

Compressed-sensing based approximation of solutions to high-dimensional operator equations

-----------------------------------
This toolbox is provided as is so far. There is absolutely no guarantees that it will work on any of your problems!
-----------------------------------

DISCLAIMER: This is a python package relying on the FEniCS solver. I will *not* provide assistance with installing FEniCS (but will happily help with your particular problem regarding this implementation, ellipticity conditions, or any compressed sensing remarks)


About this project
==================

This project stems from the collaboration between 
* Dr. Jean-Luc Bouchot, Chair C for Mathematics (Analysis), RWTH Aachen ( bouchot@mathc.rwth-aachen.de )
* Benjamin Bykowski, RWTH Aachen
* Prof. Holger Rauhut, Chair C for Mathematics (Analysis), RWTH Aachen
* Prof. Christoph Schwab, Seminar for Applied Mathematics, ETH Zurich


Requirements
============
As we are not professional software engineers, take these considerations with caution. 
Our implementation has been tested on the High Performance Cluster of the RWTH Aachen (more details on the architecture may be found here: http://www.itc.rwth-aachen.de/cms/IT_Center/IT_Center/Aktuelle_Meldungen/~fehc/Wir_heissen_jetzt_IT_Center/?lidx=1 )
To sum up, here is what you may need:
* FEniCS: http://fenicsproject.org/ ; used as the black box solver for the PDEs
* ProgressBar: https://pypi.python.org/pypi/progressbar ; used to keep track of the advancement of the PDE solves
* Scientific Linux (including in particular scipy, numpy, etc... ) - Ubuntu works just fine too!


Accompanying papers - Theory
============================
The details of these methods can be found in the following publications (pdf will follow):
* H. Rauhut and C. Schwab 
"Compressive sensing Petrov-Galerkin approximation of high-dimensional parametric operator equations"
in revision

* J.-L. Bouchot, B. Bykowski, H. Rauhut and C. Schwab
"Compressed sensing Petrov-Galerkin approximations for parametric PDEs"
Sampling Theory and Applications 2015 (SampTA 15)

* B. Bykowski
"Weigthed l1 minimization methods for numerical approximations of parametric PDEs under uncertainty quantification"
Master's Thesis, Chair C for Mathematics, RWTH Aachen, July 2015

* J.-L. Bouchot, R. Rauhut and C. Schwab
"A multi-level compressed sensing Petrov Galerkin method for the approximation of solution to high dimensional parametric operators"
In preparation (estimated submission September 2015)


License
=======
More details very soon.

Test examples
=============
Still to come.