CSPDEs
======

Compressed-sensing based approximation of solutions to high-dimensional operator equations

-----------------------------------
This toolbox is provided as is so far. There is absolutely no guarantees that it will work on any of your problems!
-----------------------------------

-----------------------------------

DISCLAIMER: This is a python package relying on the FEniCS solver. I will *not* provide assistance with installing FEniCS (but will happily help with your particular problem regarding this implementation, ellipticity conditions, or any compressed sensing remarks)


About this project
==================

This project stems from the collaboration between 
* Dr. Jean-Luc Bouchot, School of Mathematics and Statistics, Beijing Institute of Technology ( jlbouchot@bit.edu.cn )
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
* cvxpy: http://www.cvxpy.org/en/latest/install/index.html ; used for the convex minimizations needed for the compressed sensing part


Package Installation
====================

We recommend using anaconda for most of the installation. 
Be aware that these packages are evolving very fast and you might need to adjust here and there some versions of the various packages (most recently (2019/06/04): FEniCS and the latest version of cvxpy are not compatible). 
Try the following sequence of commands: 

First, create a separate environment for the project (note that you could use the [env-fenics.yml file](env-fenics.yml))
```
conda create --name env-fenics 
conda activate env-fenics
```

Make sure you have `FEniCS` installed
```
conda install -c conda-forge fenics
```

Add `cvxpy` with the right version
```
conda config --add channels oxfordcontrol
conda install -c conda-forge lapack
conda install -c cvxgrp cvxpy=1.0.11
```

(Note that the specified version is not the latest one!)
```
conda install nose
nosetests cvxpy
```

`numba` can be useful for certain speed-ups (not implemented yet)
```
conda install numba
```

You also need to install ```progressbar``` separately using -- otherwise, you might be standing in front of your computer, not knowing what is happening
```
python setup.py install 
```
in the folder where you downloaded and extracted the archive. 






### Installation
Nothing difficult here: simply clone this repository: 
```
git clone https://github.com/jlbouchot/CSPDEs.git
cd CSPDE
```
and you are ready to start playing with things. 



Accompanying papers - Theory
============================
The details of these methods can be found in the following publications:
* H. Rauhut and C. Schwab 
"Compressive sensing Petrov-Galerkin approximation of high-dimensional parametric operator equations"
Mathematics of Computation 86(304):661-700, 2017, 
([Preprint](http://www.mathc.rwth-aachen.de/~rauhut/files/csparampde.pdf))

* J.-L. Bouchot, B. Bykowski, H. Rauhut and C. Schwab
"Compressed sensing Petrov-Galerkin approximations for parametric PDEs"
Sampling Theory and Applications 2015 (SampTA 15), ([Preprint](http://www.mathc.rwth-aachen.de/~rauhut/files/SampTA15_BBRS.pdf))

* B. Bykowski
"Weigthed l1 minimization methods for numerical approximations of parametric PDEs under uncertainty quantification"
Master's Thesis, Chair C for Mathematics, RWTH Aachen, July 2015. ([Thesis](./Papers/bykowski_master.pdf))

* J.-L. Bouchot, R. Rauhut and C. Schwab
"Multi-level Compressed Sensing Petrov-Galerkin discretization of high-dimensional parametric PDEs"
Submitted Jan. 2017, ([Preprint](http://www.mathc.rwth-aachen.de/~rauhut/files/MLCSPG.pdf))

* J.-L. Bouchot
"Weighted block compressed sensing for parametrized function approximation"
Submitted Nov. 2018, ([Preprint](https://arxiv.org/abs/1811.04598))

License
=======
More details very soon.

Test examples
=============
You will find some (dirty) test examples, usually named test_XXX. In the notation, diff stands for a diffusion problem, ML for Multi-level, v for constant coefficients while p is for polynomially growing weights, pwld for PieceWise Linear Diffusion. 
Better structured tests will be made available soon. 

Acknowledgments
===============
This work was partly supported by the European Research Council through the grant StG 258926. Part of this work was developed as J.-L. B. and H. R. were visiting the Hausdorff Research Center for Mathematics as part of the Hausdorff Trimester Program on Mathematics of Signal Processing. 
