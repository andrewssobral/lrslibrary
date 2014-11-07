function sz = size(a,idx)
%SIZE Size of tenmat.
%
%   D = SIZE(X) returns the two-element row vector D = [M N]
%   containing the number of rows and columns in the matrix.
% 
%   M = SIZE(X,DIM) returns the length of the dimension specified by
%   the scalar DIM.  For example, SIZE(X,1) returns the number of
%   rows.
%
%   See also TENMAT, TENMAT/TSIZE.
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


if isempty(a.data)
    sz = [];
elseif exist('idx', 'var')
    sz = size(a.data, idx);
else
    sz = size(a.data);
end
