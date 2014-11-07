function C = mtimes(A,B)
%MTIMES sptensor-scalar multiplication.
% 
%   C = MTIMES(A,B) is called for the syntax 'A * B' when A or B is a
%   sparse tensor and the other argument is a scalar.
% 
%   For tensor-matrix multiplication, use TTM.
%   For tensor-tensor multiplication, use TTT.
%   For tensor-tensor array multiplication, use TIMES or 'A .* B'.
%
%   See also SPTENSOR, SPTENSOR/TTM, SPTENSOR/TTT, SPTENSOR/TIMES
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


if isscalar(B)    
    C = sptensor(A.subs, A.vals * B, size(A));
    return;
end

if isscalar(A)
    C = sptensor(B.subs, B.vals * A, size(B));
    return;
end

error('MTIMES only supports the scalar case for sparse tensors');
