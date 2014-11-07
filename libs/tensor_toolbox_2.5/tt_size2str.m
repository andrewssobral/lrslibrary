function s = tt_size2str(sz)
%TT_SIZE2STR Convert size to a string that can be printed.
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


if isempty(sz)
    s = sprintf('[empty tensor]');
    return;
end

if numel(sz) == 1
    s = sprintf('%d',sz);
else
    s = [sprintf('%d x ',sz(1:end-1)) sprintf('%d', sz(end)) ];
end

