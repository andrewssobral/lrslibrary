function [ Yn,hist] = lraNTD_ANLS( Y,opts )
%% Nonnegative Tucker Decomposition based on Low-rank Approximation speedup.
%   Usage: Yn = lraNTD_ANLS( Y,opts );
%  Output: Yn is a ttensor
%    opts.
%         NumOfComp: a vector specifying the dimension of Yn.core;
%         nlssolver: [HALS]|APG|MU the solver to be used for nonnegative least-squares
%         corealg: [APG]|MU|none: algorithm for nonnegative core. 'none'
%            will return a real valued core tensor instead of an nonnegative
%            one.
%         maxiter: [100] max number of iterations
%         maxiniter: [20] max number of iterations for internal loop (for each
%               sub-nls problem)
%         tol: 1e-6 the algorithm terminates if ||A(1)-A(1)_old||<tol
%         trackit: [20] check the results after each 'trackit' iterations
%         tdalgFile: a string specifying the file for unconstrained
%               Tucker decomposition. Otherwise 'tucker_als' with 5
%               iterations will be called. 
%         sparsity: an (N+1)-by-1 nonnegative vector whose entries should be
%                significatnly less than 1. sparsity(N+1) is used for the
%                sparsity of the core tensor.
%   This code depends on the TensorToolbox which is available at:
%           http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.5.html
%
%  If you think this algorithm is useful, please cite
%    Guoxu Zhou; Cichocki, A.; Shengli Xie; , "Fast Nonnegative Matrix/Tensor Factorization Based on Low-Rank Approximation,"
%    IEEE Transactions on Signal Processing, vol.60, no.6, pp.2928-2940, June 2012
%    doi: 10.1109/TSP.2012.2190410
%    URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6166354&isnumber=6198804
%
%   by Guoxu Zhou
%   http://www.bsp.brain.riken.jp/~zhougx/tensor.html
%   E-mail: zhouguoxu@gmail.com or zhouguoxu@ieee.org
%
%   Last updated: Sept-5, 2013
