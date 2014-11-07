function u = nvecs(X,n,r,opts)
%NVECS Compute the leading mode-n vectors for a ttensor.
%
%   U = NVECS(X,n,r) computes the r leading eigenvalues of Xn*Xn'
%   (where Xn is the mode-n matricization of X), which provides
%   information about the mode-n fibers. In two-dimensions, the r
%   leading mode-1 vectors are the same as the r left singular vectors
%   and the r leading mode-2 vectors are the same as the r right
%   singular vectors.
%
%   U = NVECS(X,n,r,OPTS) specifies options:
%   OPTS.eigsopts: options passed to the EIGS routine [struct('disp',0)]
%   OPTS.flipsign: make each column's largest element positive [true]
%
%   See also TTENSOR, TENMAT, EIGS.
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
% $Id: nvecs.m,v 1.7 2010/03/19 23:46:31 tgkolda Exp $

if ~exist('opts','var')
    opts = struct;
end

if isfield(opts,'eigsopts')
    eigsopts = opts.eigsopts;
else
    eigsopts.disp = 0;
end

% Compute inner product of all n-1 factors
V = cell(ndims(X),1);
for i = 1:ndims(X)
    if i == n, 
        V{i} = X.u{i};
    else
        V{i} = X.u{i}' * X.u{i};
    end
end

% Form H
H = ttm(X.core,V);

if isa(H,'sptensor')
    HnT = double(sptenmat(H,n,'t'));
else
    H = full(H);
    HnT = double(tenmat(H,n,'t'));
end
G = X.core;
if isa(G,'sptensor')
    GnT = double(sptenmat(G,n,'t'));
else
    G = full(G);
    GnT = double(tenmat(G,n,'t'));
end

% Compute Xn * Xn'
Y = HnT'*GnT*X.u{n}';

[u,d] = eigs(Y, r, 'LM', eigsopts);

if isfield(opts,'flipsign') 
    flipsign = opts.flipsign;
else
    flipsign = true;
end
    
if flipsign
    % Make the largest magnitude element be positive
    [val,loc] = max(abs(u));
    for i = 1:r
        if u(loc(i),i) < 0
            u(:,i) = u(:,i) * -1;
        end
    end
end
