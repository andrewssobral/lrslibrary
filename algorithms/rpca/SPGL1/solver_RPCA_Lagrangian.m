function [varargout] = solver_RPCA_Lagrangian(Y,lambda_L,lambda_S,varargin)
% [L,S,errHist] = solver_RPCA_Lagrangian(Y,lambda_L,lambda_S, A_cell, opts)
% Solves the problem
%   minimize_{L,S} .5|| L + S - Y ||_F^2 + lambda_L ||L||_* + lambda_S ||S||_1
%
%   or if A_cell is provided, where A_cell = {A, At}
%   (A is a function handle, At is a function handle to the transpose of A)
%   then
%
%   minimize_{L,S} .5|| A(L + S) - Y ||_F^2 + lambda_L ||L||_* + lambda_S ||S||_1
%   (here, Y usually represents A(Y) )
%
%   errHist(1,:) is a record of the residual
%   errHist(2,:) is a record fo the full objective
%   errHist(3,:) is the output of opts.errFcn if provided
%
% opts is a structure with options:
%   opts.sum, opts.max  (as described above)
%   opts.L0         initial guess for L (default is 0)
%   opts.S0         initial guess for S (default is 0)
%   opts.tol        sets stopping tolerance
%   opts.maxIts     sets maximum number of iterations
%   opts.printEvery will print information this many iterations
%   opts.displayTime will print out timing information (default is true for large problems)
%   opts.errFcn     a function of (L,S) that records information
%   opts.trueObj    if provided, this will be subtracted from errHist(2,:)
%   opts.Lip        Lipschitz constant, i.e., 2*spectralNorm(A)^2
%                       by default, assume 2 (e.g., good if A = P_Omega)
%   opts.FISTA      whether to use FISTA or not. By default, true
%     opts.restart  how often to restart FISTA; set to -Inf to make it automatic
%   opts.BB         whether to use the Barzilai-Borwein spectral steplength
%     opts.BB_type  which BB stepsize to take. Default is 1, the larger step
%     opts.BB_split whether to calculate stepslengths for S and L independently.
%       Default is false, which is recommended.
%   opts.quasiNewton  uses quasi-Newton-like Gauss-Seidel scheme.
%                     Only available in "max" mode
%     opts.quasiNewton_stepsize     stepsize length. Default is .8*(2/Lip)
%     opts.quasinewton_SLS          whether to take S-L-S sequence (default is true)
%                                   otherwise, takes a L-S Gauss-Seidel sequence
%   opts.SVDstyle   controls what type of SVD is performed.
%       1 = full svd using matlab's "svd". Best for small problems
%       2 = partial svd using matlab's "svds". Not recommended.
%       3 = partial svd using PROPACK, if installed. Better than option 2, worse than 4
%       4 = partial svd using randomized linear algebra, following
%           the Halko/Tropp/Martinnson "Structure in Randomness" paper
%       in option 4, there are additional options:
%       opts.SVDwarmstart   whether to "warm-start" the algorithm
%       opts.SVDnPower  number of power iterations (default is 2 unless warm start)
%       opts.SVDoffset  oversampling, e.g., "rho" in Tropp's paper. Default is 5
%
% Stephen Becker, March 6 2014
% See also solver_RPCA_constrained.m


[varargout{1:nargout}] = solver_RPCA_constrained(Y,lambda_S/lambda_L, -lambda_L, varargin{:});