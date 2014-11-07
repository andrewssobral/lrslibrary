function t = ktensor(varargin)
%KTENSOR Tensor stored as a Kruskal operator (decomposed).
%
%   K = KTENSOR(lambda,U1,U2,...,UM) creates a Kruskal tensor from its
%   constituent parts. Here lambda is a k-vector and each Um is a
%   matrix with k columns.
%
%   K = KTENSOR(lambda, U) is the same as above except that U is a
%   cell array containing matrix Um in cell m.
%
%   K = KTENSOR(U) assumes U is a cell array containing matrix Um in
%   cell m and assigns the weight of each factor to be one.
%
%   K = KTENSOR(T) creates a ktensor by copying an existing ktensor.
%
%   Examples
%   K = ktensor([3; 2], rand(4,2), rand(5,2), rand(3,2))
%
%   See also TENSOR, TTENSOR, KTENSOR/FULL.
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


% EMPTY CONSTRUCTOR
if nargin == 0
    t.lambda = [];
    t.u = {};
    t = class(t,'ktensor');
    return;
end

% Copy CONSTRUCTOR
if (nargin == 1) && isa(varargin{1}, 'ktensor')
    t.lambda = varargin{1}.lambda;
    t.u = varargin{1}.u;
    t = class(t, 'ktensor');
    return;
end

if isa(varargin{1},'cell')

    u = varargin{1};
    t.lambda = ones(size(u{1},2),1);
    t.u = u;
    
else

    t.lambda = varargin{1};
    if ~isa(t.lambda,'numeric') || ndims(t.lambda) ~=2 || size(t.lambda,2) ~= 1
	error('LAMBDA must be a column vector.');
    end
    
    if isa(varargin{2},'cell')
	t.u = varargin{2};
    else
	for i = 2 : nargin
	    t.u{i-1} = varargin{i};
	end
    end

end
    
    
% Check that each Um is indeed a matrix
for i = 1 : length(t.u)
    if ndims(t.u{i}) ~= 2
	error(['Matrix U' int2str(i) ' is not a matrix!']);
    end
end

% Size error checking			     
k = length(t.lambda); 
for i = 1 : length(t.u)            
    if  size(t.u{i},2) ~= k
       error(['Matrix U' int2str(i) ' does not have ' int2str(k) ' columns.']);
    end
end

t = class(t, 'ktensor');
return;
