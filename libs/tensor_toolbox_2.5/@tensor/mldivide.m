function Z = mldivide(X,Y)
%MLDIVIDE Slash left division for tensors.
%
%   MLDIVIDE(A,B) is called for the syntax 'A \ B' when A is a scalar and B
%   is a tensor.  
%
%   Example
%   X = tenrand([4 3 2],5);
%   3 \ X
%
%   See also TENSOR, TENSOR/LDIVIDE.
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


if isscalar(X)
    Z = tenfun(@ldivide,X,Y);
    return;
end

error('MLDIVIDE only supports the scalar case for tensors');

