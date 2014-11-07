function C = times(A,B)
%TIMES Element-wise multiplication for ktensor.
%
%   TIMES(A,B) denotes element-by-element multiplication.  
% 
%   C = TIMES(A,B) is called for the syntax 'A .* B' when A or B is a
%   tensor.
%
%   See also KTENSOR, SPTENSOR/TIMES.
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


if ~isequal(size(A),size(B))
    error('Must be two tensors of the same size');
end

switch class(B)
    case {'sptensor','tensor'}
        % Call back to sptensor version.
        C = times(B,A);
        return;
    otherwise
        error('Invalid second argument for ktensor/times');
end
