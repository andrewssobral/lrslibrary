function C = mtimes(A,B)
%MTIMES Multiplies two tenmat objects.
%
%  C = MTIMES(A,B) computes the product of A and B. The result is a
%  TENMAT object and can be transformed into a tensor.
%
%  C = MTIMES(A,B) is called for the syntax 'A * B' when A or B is a
%  TENMAT object.  
%
%   See also TENMAT.
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


% Handle scalar input
if ~isa(B,'tenmat') && numel(B) == 1
    C = A;
    C.data = C.data * B;
    return;
end
if ~isa(A,'tenmat') && numel(A) == 1
    C = B;
    C.data = C.data * A;
    return;
end

% Handle matrix input
if ~isa(A,'tenmat')
    A = tenmat(A,1);
end

if ~isa(B,'tenmat')
    B = tenmat(B,1);
end

% Error check
if size(A,2) ~= size(B,1)  
    error(['Size mismatch: Number of columns in A is not equal to' ...
	   ' the number of rows in B']);
end

tsiz = [A.tsize(A.rindices) B.tsize(B.cindices)];

if ~isempty(tsiz)
    C = tenmat;
    C.tsize = tsiz;
    C.rindices = 1:length(A.rindices);
    C.cindices = (1:length(B.cindices)) + length(A.rindices);
    C.data = A.data * B.data;
else
    C = A.data * B.data;
end

