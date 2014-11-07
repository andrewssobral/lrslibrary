function ok = tt_subscheck(subs)
%TT_SUBSCHECK Checks for valid subscripts.
%
%  TT_SUBSCHECK(S) throws an error if S is not a valid subscript
%  array, which means that S is a matrix of real-valued, finite,
%  positive, integer subscripts.
%
%  X = TT_SUBSCHECK(S) returns true if S is a valid and false
%  otherwise.
%
%  See also TT_SIZECHECK, TT_VALSCHECK.
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

%
% Includes improvements offered by Marcus Brubaker.

if isempty(subs)
    ok = true;
elseif ndims(subs) == 2 && isreal(subs) ...
        && all(isfinite(subs(:)) & subs(:) > 0) ...
        && isequal(subs,round(subs)) 
    ok = true;
else
    ok = false;
end

if ~ok && nargout == 0
    error('Subscripts must be a matrix of real positive integers');
end
