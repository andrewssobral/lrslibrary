function ok = tt_sizecheck(siz)
%TT_SIZECHECK Checks that the size is valid.
%
%  TT_SIZECHECK(S) throws an error if S is not a valid size array,
%  which means that it is a row vector with strictly postitive,
%  real-valued, finite integer values.
%
%  X = TT_SIZECHECK(S) returns true if S is a valid and false otherwise.
%
%  See also TT_SUBSCHECK.
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


if ndims(siz) == 2 && size(siz,1) == 1 ...
        && isreal(siz) ...
        && ~any(isnan(siz(:))) && ~any(isinf(siz(:))) ...
        && isequal(siz,round(siz)) && all(siz(:) > 0) 
    ok = true;
else
    ok = false;
end

if ~ok && nargout == 0
    error('Size must be a row vector of real positive integers');
end
