function res = innerprod(X,Y)
%INNERPROD Efficient inner product with a ktensor.
%
%   R = INNERPROD(X,Y) efficiently computes the inner product between
%   two tensors X and Y.  If Y is a ktensor, the inner product is
%   computed using inner products of the factor matrices, X{i}'*Y{i}.
%   Otherwise, the inner product is computed using ttv with all of
%   the columns of X's factor matrices, X{i}.
%
%   See also KTENSOR, KTENSOR/TTV
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


if ~isequal(size(X),size(Y))
    error('X and Y must be the same size.');
end

% X is a ktensor
switch class(Y)
 
  case {'ktensor'}
    M = X.lambda * Y.lambda';
    for n = 1:ndims(X)
        M = M .* (X.u{n}' * Y.u{n});
    end
    res = sum(M(:));
    
  case {'tensor','sptensor','ttensor'}
    R = length(X.lambda);
    vecs = cell(1,ndims(X));
    res = 0;
    for r = 1:R
      for n = 1:ndims(X)
        vecs{n} = X.u{n}(:,r);
      end
      res = res + X.lambda(r) * ttv(Y,vecs);
    end
    
  otherwise
    disp(['Inner product not available for class ' class(Y)]);
end

return;
