function t = ttensor(varargin)
%TTENSOR Tensor stored as a Tucker operator (decomposed).
%
%   T = TTENSOR(G,U1,U2,...,UM) creates a TUCKER tensor from its
%   constituent parts. Here G is a tensor of size K1 x K2 x ... x KM
%   and each Um is a matrix with Km columns.
%
%   T = TTENSOR(G,U) is the same as above except that U is a cell
%   array containing matrix Um in cell m.
%
%   The core tensor G can be any type of tensor that supports the
%   following functions:
%     - size
%     - uminus
%     - disp (with 2 arguments; see, e.g., TENSOR/DISP)
%     - ttv
%     - ttm
%     - mtimes (scalar multiplication only)
%     - permute
%     - subsasgn
%     - subsref
%   
%   T = TTENSOR(S) creates a TUCKER tensor by copying an existing
%   TUCKER tensor.
%
%   T = TTENSOR is the empty constructor.
%
%   See also TENSOR, KTENSOR
%
%MATLAB Tensor Toolbox.
%Copyright 2010, Sandia Corporation. 

% This is the MATLAB Tensor Toolbox by Brett Bader and Tamara Kolda. 
% http://csmr.ca.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2010) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in tensor_toolbox/LICENSE.txt
% $Id: ttensor.m,v 1.7 2010/03/19 23:46:31 tgkolda Exp $

% Empty constructor
if (nargin == 0)
    t.core = tensor;                    % empty tensor
    t.u = [];
    t = class(t, 'ttensor');
    return;
end

% Copy CONSTRUCTOR
if (nargin == 1) && isa(varargin{1}, 'ttensor')
    t.core = varargin{1}.core;
    t.u = varargin{1}.u;
    t = class(t, 'ttensor');
    return;
end

% Core can be basically anything that supports certain functions.
t.core = varargin{1};

if isa(varargin{2},'cell')
    t.u = varargin{2};
else
    for i = 2 : nargin
	t.u{i-1} = varargin{i};
    end
end

% Check that each Um is indeed a matrix
for i = 1 : length(t.u)
    if ndims(t.u{i}) ~= 2
	error(['Matrix U' int2str(i) ' is not a matrix!']);
    end
end

% Size error checking			     
k = size(t.core); 

if length(k) ~= length(t.u)
    error(['CORE has order ', int2str(length(k)), ...
	   ' but there are ', int2str(length(t.u)), ' matrices.']);
end

for i = 1 : length(t.u)            
    if  size(t.u{i},2) ~= k(i)
	error(['Matrix U' int2str(i) ' does not have ' int2str(k(i)) ...
	       ' columns.']);
    end
end

t = class(t, 'ttensor');
return;
