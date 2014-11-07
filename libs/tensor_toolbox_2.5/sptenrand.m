function Y = sptenrand(sz,nz)
%SPTENRAND Sparse uniformly distributed random tensor.
%
%   R = SPTENRAND(sz,density) creates a random sparse tensor of the
%   specified sz with approximately density*prod(sz) nonzero
%   entries.
%
%   R = SPTENRAND(sz,nz) creates a random sparse tensor of the
%   specified sz with approximately nz nonzero entries.
%
%   Example: R = sptenrand([5 4 2],12);
%
%   See also SPTENSOR, TENRAND, RAND.
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



% Error check on siz
if ndims(sz) ~= 2 || size(sz,1) ~= 1
    error('Size must be a row vector');
end

% Error check on nz
if ~exist('nz','var') || (nz < 0)
    error('2nd argument must be positive');
end

% Is nz an number or a fraction? Ultimately want a number.
if (nz < 1)
    nz = prod(sz) * nz;
end

% Make sure nz is an integer
nz = ceil(nz);

% Keep iterating until we find enough unique nonzeros or we
% give up
subs = [];
cnt = 0;
while (size(subs,1) < nz) && (cnt < 10)
    subs = ceil( rand(nz, size(sz,2)) * diag(sz) );
    subs = unique(subs, 'rows');
    cnt = cnt + 1;
end

% Extract nnz subscipts and create a corresponding list of
% values
nz = min(nz, size(subs,1));
subs = subs(1:nz,:);
vals = rand(nz,1);

Y = sptensor(subs,vals,sz);
return;
