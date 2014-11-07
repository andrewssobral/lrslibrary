function C = mtimes(A,B)
%MTIMES tensor-scalar multiplication.
% 
%   C = MTIMES(A,B) is called for the syntax 'A * B' when A or B is a
%   tensor and the other argument is a scalar.
% 
%   For tensor-matrix multiplication, use TTM.
%   For tensor-tensor multiplication, use TTT.
%   For tensor-tensor array multiplication, use TIMES or 'A .* B'.
% 
%   Examples
%   X = tenrand([3,4,2])
%   W = 5 * X
%
%   See also TENSOR, TENSOR/TTM, TENSOR/TTT, TENSOR/TIMES
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


%%
if isscalar(B)
    C = A;
    C.data = B * C.data;
    return;
end

if isscalar(A)
    C = B;
    C.data = A * C.data;
    return;
end

error('Mtimes only supports a tensor times a scalar');





