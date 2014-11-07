function f = tt_loglikelihood(X,M)
%TT_LOGLIKELIHOOD Compute log-likelihood of data X with model M.
%
%   F = TT_LOGLIKELIHOOD(X,M) computes the log-likelihood of model M given
%   data X, where M is a ktensor and X is a tensor or sptensor.
%   Specifically, F = - (sum_i m_i - x_i * log_i) where i is a multiindex
%   across all tensor dimensions.
%
%   See also cp_apr, tensor, sptensor, ktensor.
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

N = ndims(X);

if ~isa(M, 'ktensor')
    error('M must be a ktensor');
end

M = normalize(M,1,1);

if isa(X, 'sptensor')
    xsubs = X.subs;
    A = M.U{1}(xsubs(:,1),:);
    for n = 2:N
       A = A .* M.U{n}(xsubs(:,n),:); 
    end
    f = sum(X.vals .* log(sum(A,2))) - sum(sum(M.U{1}));
else
    f = sum(sum(double(tenmat(X,1)) .* log(double(tenmat(M,1))))) - sum(sum(M.U{1}));
end


