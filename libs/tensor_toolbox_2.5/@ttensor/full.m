function X = full(T)
%FULL Convert a ttensor to a (dense) tensor.
%
%   X = FULL(T) converts ttensor T to (dense) tensor X.
%
%   See also TTENSOR, TENSOR.
%
%MATLAB Tensor Toolbox.
%Copyright 2010, Sandia Corporation. 

% This is the MATLAB Tensor Toolbox by Brett Bader and Tamara Kolda. 
% http://csmr.ca.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2010) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in tensor_toolbox/LICENSE.txt
% $Id: full.m,v 1.11 2010/03/19 23:46:31 tgkolda Exp $

% Preallocate to ensure there is enough space
X = tenzeros(size(T));

% Now do the calculation 
X = ttm(T.core,T.u);

% Make sure that X is a dense tensor (small chance it could be a sparse
% tensor).
X = tensor(X);

return;
