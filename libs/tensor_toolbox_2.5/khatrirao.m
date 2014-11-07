function P = khatrirao(varargin)
%KHATRIRAO Khatri-Rao product of matrices.
%
%   KHATRIRAO(A,B) computes the Khatri-Rao product of matrices A and
%   B that have the same number of columns.  The result is the
%   column-wise Kronecker product
%   [KRON(A(:,1),B(:,1)) ... KRON(A(:,n),B(:,n))]
%
%   KHATRIRAO(A1,A2,...) computes the Khatri-Rao product of
%   multiple matrices that have the same number of columns.
%
%   KHATRIRAO(C) computes the Khatri-Rao product of
%   the matrices in cell array C.
%
%   KHATRIRAO(...,'r') computes the Khatri-Rao product in reverse
%   order.
%
%   Examples
%   A = rand(5,2); B = rand(3,2); C = rand(2,2);
%   khatrirao(A,B) %<-- Khatri-Rao of A and B
%   khatrirao(B,A,'r') %<-- same thing as above
%   khatrirao({C,B,A}) %<-- passing a cell array
%   khatrirao({A,B,C},'r') %<-- same as above
%
%   See also TENSOR, KTENSOR.
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


%% Error checking on input and set matrix order
% Note that this next if/else check forces A to be a cell array.
if ischar(varargin{end}) && varargin{end} == 'r'
    if nargin == 2 && iscell(varargin{1})
        % Input is a single cell array
        A = varargin{1};
    else
        % Input is a sequence of matrices
        A = {varargin{1:end-1}};
    end
    matorder = length(A):-1:1;
else
    if nargin == 1 && iscell(varargin{1})
        % Input is a single cell array
        A = varargin{1};
    else
        % Input is a sequence of matrices
        A = varargin;
    end
    matorder = 1:length(A);
end

%% Error check on matrices and compute number of rows in result 

% N = number of columns (must be the same for every input)
N = size(A{1},2); 

% After loop, M = number of rows in the result
M = 1; 

for i = matorder
    if ndims(A) ~= 2
        error('Each argument must be a matrix');
    end
    if (N ~= size(A{i},2))
        error('All matrices must have the same number of columns.')
    end
    M = M * size(A{i},1);
end

%% Computation

% Preallocate
P = zeros(M,N);

% Loop through all the columns
for n = 1:N
    % Loop through all the matrices
    ab = A{matorder(1)}(:,n);
    for i = matorder(2:end)
       % Compute outer product of nth columns
       ab = A{i}(:,n) * ab(:).';
    end
    % Fill nth column of P with reshaped result
    P(:,n) = ab(:);          
end

