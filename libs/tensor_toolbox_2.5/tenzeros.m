function X = tenzeros(varargin)
%TENZEROS Create zeros tensor.
%
%   X = TENZEROS(SZ) forms a tensor of size SZ with all zeros.
%
%   TENZEROS(SZ) is equivalent to TENSOR(ZEROS(SZ(1),SZ(2),...),SZ).
%
%   See also TENSOR, ZEROS.
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


if nargin == 1
    sz = varargin{1};
else
    sz = cell2mat(varargin);
end

if isempty(sz)
    X = tensor;
    return;
end

if nargin == 2
    order = sz;
    dim = varargin{1};
    sz = dim * ones(1,order);
end

data = zeros([sz 1 1]);
X = tensor(data,sz);

