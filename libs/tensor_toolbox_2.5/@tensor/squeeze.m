function Y = squeeze(X)
%SQUEEZE Remove singleton dimensions from a tensor.
%
%   Y = SQUEEZE(X) returns a tensor Y with the same elements as
%   X but with all the singleton dimensions removed.  A singleton
%   is a dimension such that size(X,dim)==1.  
%
%   If X has *only* singleton dimensions, then Y is a scalar.
%
%   Examples
%   squeeze( tenrand([2,1,3]) ) %<-- returns a 2-by-3 tensor
%   squeeze( tenrand([1 1]) ) %<-- returns a scalar
%
%   See also TENSOR.
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


if all(X.size > 1)
  % No singleton dimensions to squeeze
  Y = X;
else
  idx = find(X.size > 1);
  if numel(idx) == 0
    % Scalar case - only singleton dimensions
    Y = X.data;
  else
    siz = X.size(idx);
    Y = tensor(squeeze(X.data),siz);
  end
end

return;
