function t = ones(t)
%ONES Replace nonzero elements of sparse tensor with ones.
%
%   S = ONES(T) generates a sparse tensor with the same sparsity
%   structure as T, but with ones in the nonzero position.
%
%   See also SPTENSOR, SPONES.
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


t.vals = ones(size(t.vals));
