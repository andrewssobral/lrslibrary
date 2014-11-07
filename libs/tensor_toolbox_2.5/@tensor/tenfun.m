function Z = tenfun(fun,varargin)
%TENFUN Apply a function to each element in a tensor.
%
%   TENFUN(F,X,...) applies the function specified by the function
%   handle F to the given arguments.  Either both arguments
%   must be tensors, or one is a tensor and the other is a scalar/MDA.
%
%   Examples
%   Z = tenfun(@(x)(x+1),X) %<-- increase every element by one
%   Z = tenfun(@eq,X,1) %<-- logical comparison of X with scalar
%   Z = tenfun(@plus,X,Y) %<-- adds the two tensors X and Y.
%   Z = tenfun(@max,X,Y,Z) %<-- max over all elements in X,Y,Z
%
%   See also TENSOR.
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


if nargin < 2
    %  error('TENFUN requires at least two input arguments')
    error('Not enough input arguments.');
end

if ~isa(fun, 'function_handle')
    error('First argument must be a function handle.');
end

%% Case I: NARGIN == 2 (one function and one tensor)
if (nargin == 2) && isa(varargin{1},'tensor')
    Z = varargin{1};
    Z.data(:) = fun(Z.data(:));
    return;
end

%% Determine if function is binary 
% Note that we swap the arguments below if 2nd argument is sparse, but need
% to take special measures for those functions that don't commute.
binfuns = {@plus,@minus,@eq,@ge,@gt,@le,@lt,@ne,@and,@or,@xor, ...
           @power,@times,@ldivide,@rdivide};
isbinary = false;
for i = 1 : numel(binfuns)
    if isequal(fun,binfuns{i})
        isbinary = true;
        break;
    end
end

%% Case II: NARGIN == 3 and function is binary
if (nargin == 3) && isbinary
    
    X = varargin{1};
    Y = varargin{2};
    
    % Case IIa: X is a scalar/MDA and Y is a tensor
    if isnumeric(X) && isa(Y,'tensor')
        Z = Y;
        Z.data = fun(X,Y.data);
        return;
    end

    % Case IIb: X is a tensor and Y is a scalar/MDA
    if isa(X,'tensor') && isnumeric(Y)
        Z = X;
        Z.data = fun(X.data,Y);
        return;
    end

    % Case IIc: X and Y are both tensors
    if isa(X,'tensor') && isa(Y,'tensor')
        if ~(isequal(size(X),size(Y)))
            error('Tensor size mismatch.')
        end
        data = fun(X.data,Y.data);
        Z = tensor(data,size(X));
        return;
    end
    
    % Case IId: Either X or Y is a sptensor
    if isa(X,'tensor') && isa(Y,'sptensor')
        if isequal(fun,@lt)
            Z = gt(Y,X);
        elseif isequal(fun,@le)
            Z = ge(Y,X);
        elseif isequal(fun,@gt)
            Z = lt(Y,X);
        elseif isequal(fun,@ge)
            Z = le(Y,X);
        elseif isequal(fun,@minus)
            Z = plus(-Y,X);
        elseif isequal(fun,@power)
            error('Cannot do array power with a sparse and dense tensor');
        elseif isequal(fun,@ldivide)
            error('Cannot do ldivide with a sparse and dense tensor');
        elseif isequal(fun,@rdivide)
            error('Cannot do rdivide with a sparse and dense tensor');
        else
            Z = fun(Y,X);
        end
        return;
    end
    
    % Case IIe: Either X or Y is not a tensor   
    error('For a binary function, arguments must be either two tensors or one tensor and a scalar/MDA.');
end

%% Case III: one function and the rest are same-sized tensors

X = varargin;
n = numel(X);
sz = size(X{1});
m = prod(sz);
Y = zeros(m,n);
for j = 1:n
    Y(:,j) = X{j}.data(:);
end
data = zeros(m,1);
for i = 1:m
    data(i) = fun(Y(i,:));
end
Z = tensor(data,sz);
return;


