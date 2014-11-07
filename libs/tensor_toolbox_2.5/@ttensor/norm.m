function nrm = norm(X)
%NORM Norm of a ttensor.
%
%   See also TTENSOR.
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
% $Id: norm.m,v 1.9 2010/03/19 23:46:31 tgkolda Exp $


if prod(size(X)) > prod(size(X.core))
    V = cell(ndims(X),1);
    for n = 1:ndims(X)
        V{n} = X.u{n}'*X.u{n};
    end
    Y = ttm(X.core,V);
    tmp = innerprod(Y, X.core);
    nrm = sqrt(tmp);
else
    nrm = norm(full(X));
end
