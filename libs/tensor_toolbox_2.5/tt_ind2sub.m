function subs = tt_ind2sub(siz,idx)
%TT_IND2SUB Multiple subscripts from linear indices.
%
%   SUBS = TT_IND2SUB(SIZ,INDS) returns that subscripts equivalent
%   to the linear indices in INDS for a tensor of size SIZ.
%
%   See also TT_SUB2IND, IND2SUB.
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


if isempty(idx)
    subs = [];
    return;
end

k = [1 cumprod(siz(1:end-1))];
n = length(siz);
idx = idx - 1;
for i = n : -1 : 1
    subs(:,i) = floor(idx / k(i)) + 1;
    idx = rem(idx,k(i));
end
