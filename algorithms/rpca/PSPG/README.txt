******************************************************************************************
IMPORTANT:
------------------------------------------------------------------------------------------
PROPACK package and its installer is included in the PROPACK_SVT folder.
This installer was written by Emmanuel Candès and Stephen Becker (http://www-stat.stanford.edu/~candes/svt/code.html)
Note that PROPACK_SVT folder contains .mex files that have been compiled into binaries for most common architectures/OS. 
If the currently compiled .mex files do not run on your architecture/OS properly, then you should run install_mex.m to
compile .mex file for your system.
******************************************************************************************

PSPG calling format:

[X,S,out]=pspg(D,stdev,tol,optimal_X,optimal_S)

INPUT
------------------------------------------------------------------------------------------
D: Data Matrix, where D=X*+S*+E*, X* is low-rank, S* is sparse, E* random noise
stdev: standard deivation of the components of random noise matrix
tol: stopping tolerance, see pspg.m for details
optimal_X: X* (if known)
optimal_S: S* (if known)
------------------------------------------------------------------------------------------

This folder includes following test file.

demo.m : Creates and solves randomly generated test problem