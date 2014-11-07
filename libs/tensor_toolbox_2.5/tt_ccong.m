function X = tt_ccong(m,n,gamma)
%TT_CCONG Create a random matrix with a fixed congruence
%
%   X = TT_CCONG(M,N,GAMMA) creates a matrix X of size M x N such
%   that each column of X has norm 1 and any two columns of X have an inner
%   product equal to GAMMA.
%
%   Based on code from Evrim Acar and the paper G. Tomasi and R. Bro, A
%   comparison of algorithms for fitting the PARAFAC model, Computational
%   Statistics & Data Analysis, 50: 1700-1734, 2006.
%
%MATLAB Tensor Toolbox.
%Copyright 2012, Sandia Corporation.

% This is the MATLAB Tensor Toolbox by T. Kolda, B. Bader, and others.
% http://www.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2012) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in the file LICENSE.txt

CG = gamma * ones(n,n) + (1-gamma) * eye(n);
CGR = chol(CG);
X = randn(m,n);
[Q,~] = qr(X,0);
X = Q * CGR;