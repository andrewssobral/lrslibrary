******************************************************************************************
IMPORTANT:
------------------------------------------------------------------------------------------
PROPACK package is required. The package is available at http://soi.stanford.edu/~rmunk/PROPACK/

For the convenience of the user we have included the PROPACK package and its installer in the PROPACK_SVT folder.
This installer was written by Emmanuel Candès and Stephen Becker (http://www-stat.stanford.edu/~candes/svt/code.html)
Note that PROPACK_SVT folder contains .mex files that have been compiled into binaries for most common architectures/OS. 
If the currently compiled .mex files do not run on your architecture/OS properly, then you should run install_mex.m to
compile .mex file for your system.
******************************************************************************************

NSA calling format:

[X,S,out]=nsa(D,stdev,tol,denoise_flag,optimal_X,optimal_S)

INPUT
------------------------------------------------------------------------------------------
D: Data Matrix, where D=X*+S*+E*, X* is low-rank, S* is sparse, E* random noise
stdev: standard deivation of the components of random noise matrix
tol: stopping tolerance, see nsa.m for details
denoise_flag: if denoising is on, then problem (4.6) in the paper is solved using the NSA output (X,S).
optimal_X: X* (if known)
optimal_S: S* (if known)
------------------------------------------------------------------------------------------

This folder includes following test files.

demo_1.m : Solves the problem of extracting the foreground from a noisy video
demo_2.m : Creates and solves randomly generated test problem