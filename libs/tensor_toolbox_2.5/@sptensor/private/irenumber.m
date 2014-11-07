function newsubs = irenumber(t, sz, range)
%RENUMBER indices for sptensor subsasgn
%
%  See also SPTENSOR/SUBSASGN
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


nz = nnz(t);
if (nz == 0)
    newsubs = [];
    return;
end

newsubs = t.subs;
for i = 1 : numel(range)
    r = range{i};
    if ischar(r) && r == ':'
        continue;
    elseif numel(r) == 1
        newsubs = [newsubs(:,1:i-1), r*ones(nz,1), newsubs(:,i:end)];
    else
        newsubs(:,i) = r(newsubs(:,i));
    end
end

