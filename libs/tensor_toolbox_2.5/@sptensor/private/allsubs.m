function s = allsubs(x)
%ALLSUBS Generate all possible subscripts for a sparse tensor X.
%
%   See also SPTENSOR.
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


%% Generate all possible indicies

% Preallocate (discover any memory issues here!)
s = zeros(prod(x.size),ndims(x));

% Generate appropriately sized ones vectors.
o = cell(ndims(x),1);
for n = 1:ndims(x)
    o{n} = ones(size(x,n),1);
end

% Generate each column of the subscripts in turn
for n = 1:ndims(x)
    i = o;
    i{n} = (1:size(x,n))';
    s(:,n) = khatrirao(i); 
end
