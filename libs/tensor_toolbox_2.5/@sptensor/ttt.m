function c = ttt(varargin)
%TTT Sparse tensor times sparse tensor.
% 
%   Z = TTT(X,Y) computes the outer product of tensors X and Y.
%
%   Z = TTT(X,Y,XDIMS,YDIMS) computes the contracted product of
%   tensors X and Y in the dimensions specified by the row vectors
%   XDIMS and YDIMS.  The sizes of the dimensions specified by XDIMS
%   and YDIMS must match; that is, size(X,XDIMS) must equal
%   size(Y,YDIMS).
%
%   Z = TTT(X,Y,DIMS) computes the inner product of tensors X and Y in
%   the dimensions specified by the vector DIMS.  The sizes of the
%   dimensions specified by DIMS must match; that is, size(X,DIMS)
%   must equal size(Y,DIMS).
%
%   In all cases, the result Z is a sparse tensor if it has 50% or
%   fewer nonzeros; otherwise ther result is returned as a dense
%   tensor.
%
%   Examples
%   X = sptenrand([4 2 3], 10);
%   Y = sptenrand([3 4 2], 10);
%   Z = ttt(X,Y) %<-- outer product of X and Y
%   Z = ttt(X,X,1:3) %<-- inner product of X with itself
%   Z = ttt(X,Y,[1 2 3],[2 3 1]) %<-- inner product of X & Y
%   Z = ttt(X,Y,[1 3],[2 1]) %<-- product of X & Y along specified dims
%
%   See also SPTENSOR, TENSOR/TTT, SPTENSOR/TTV, SPTENSOR/TTM.
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



%%%%%%%%%%%%%%%%%%%%%%
%%% ERROR CHECKING %%%
%%%%%%%%%%%%%%%%%%%%%%

% Check the number of arguments
if (nargin < 2)
    error('TTT requires at least two arguments.');
end

% Check the first argument
if ~isa(varargin{1}, 'sptensor')
    error('First argument must be a sptensor.');
else
    a = varargin{1};
end

% Check the second argument
if ~isa(varargin{2}, 'sptensor')
    error('Second argument must be a sptensor.');
else
    b = varargin{2};
end

% Optional 3rd argument
if nargin >= 3
    adims = varargin{3};
else
    adims = [];
end

% Optional 4th argument
if nargin >= 4
    bdims = varargin{4};
else
    bdims = adims;
end

if ~isempty(adims)
    tt_dimscheck(adims,ndims(a));
end
if ~isempty(bdims)
    tt_dimscheck(bdims,ndims(b));
end

asiz = size(a);
bsiz = size(b);
if ~isequal(asiz(adims),bsiz(bdims))
    error('Specified dimensions do not match.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPUTE THE PRODUCT %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Remaining dimensions
aremdims = setdiff(1:ndims(a),adims);
bremdims = setdiff(1:ndims(b),bdims);

if (nnz(a) == 0) || (nnz(b) == 0)
    if isempty(aremdims) && isempty(bremdims)
        c = 0;
    else
        c = sptensor([],[],[a.size(aremdims) b.size(bremdims)]);
    end
    return;
end

if isempty(adims) && isempty(bdims)
    aii = ones(nnz(a),1);
    bii = ones(nnz(b),1);
    m = 1;
else
    innersubs = [a.subs(:,adims); b.subs(:,bdims)];
    [junk1,junk2,loc] = unique(innersubs,'rows');
    aii = loc(1:nnz(a));
    bii = loc(nnz(a)+1:end);
    m = length(junk2);
end

if isempty(aremdims)
    ajj = ones(nnz(a),1);
    asubs = [];
    n = 1;
else
    [asubs,junk,ajj] = unique(a.subs(:,aremdims),'rows');
    n = length(junk);
end
if isempty(bremdims)
    bjj = ones(nnz(b),1);
    bsubs = [];
    p = 1;
else
    [bsubs,junk,bjj] = unique(b.subs(:,bremdims),'rows');
    p = length(junk);
end

aa = sparse(aii,ajj,a.vals,m,n);
bb = sparse(bii,bjj,b.vals,m,p);
cc = aa'*bb;

% Check for a scalar result, corresponding to an inner product.
if isempty(aremdims) && isempty(bremdims)
    c = sum(nonzeros(cc));
    return;
end

% If cc is a row vector, then transpose to work as a column vector
% (otherwise 'find' returns row vectors)
if size(cc,1) == 1
    [jj,ii,newvals] = find(cc');
else
    [ii,jj,newvals] = find(cc);
end

if isempty(asubs) && ~isempty(bsubs)
    newsubs = bsubs(jj,:);
elseif ~isempty(asubs) && isempty(bsubs)
    newsubs = asubs(ii,:);
else
    newsubs = [asubs(ii,:), bsubs(jj,:)];
end

c = sptensor(newsubs,newvals,[a.size(aremdims) b.size(bremdims)]);

% Convert the result to dense if it has more than 50% nonzeros.
if nnz(c) > 0.5 * prod(c.size)
    c = tensor(c);
end
