function [f,G] = tt_cp_wfg_sparse(Zvals,W,A,normZsqr,memflag)
%TT_CP_WFG_SPARSE Computes weighted CP function and gradient.
%
%   [F,G] = TT_CP_WFG_SPARSE(ZVALS,W,A) computes the function and
%   gradient with respect to A of || W .* (Z - ktensor(A)) ||^2 where
%   Z = W .* X. The variable ZVALS contains the values of the tensor Z
%   at the locations specified by W.subs. (ZVALS can be computed using
%   a provided preprocessing function.) The variable A is a cell array
%   of component matrices. The tensor W is a sparse tensor that has
%   ones in entries where we know the values.
%
%   [F,G] = TT_CP_WFG_SPARSE(ZVALS,W,A,NORMZSQR) also passes in the
%   pre-computed norm of Z, which makes the computations faster.
%
%   [F,G] = TT_CP_WFG_SPARSE(ZVALS,A,W,NORMZSQR,false) uses less memory
%   but more time and is appropriate for very large sparse tensors.
% 
%   See also TT_CP_WFG_SPARSE_SETUP, CP_WFG, CP_WFUN.
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


%% Set-up
N = ndims(W);
R = size(A{1},2);
sz = cellfun(@(x)size(x,1),A);
Wsubs = W.subs;
Wvals = W.vals;
Nvals = length(Wvals);

if ~exist('memflag','var')
    memflag = true;
end

%% Compute B = W.*ktensor(A)
Bvals = zeros(Nvals,1);
for r = 1:R
    newvals = Wvals;
    for n = 1:N
        bigArn = A{n}(Wsubs(:,n),r);
        newvals = newvals .* bigArn;
    end
    Bvals = Bvals + newvals;
end

%% Compute normZ
if ~exist('normZsqr','var')
    normZsqr = sum(Zvals.^2);
end

%% function value: f = 0.5 * normZsqr - innerprod(Z,B) + 0.5 * norm(B)^2
f = 0.5 * normZsqr - Zvals'*Bvals + 0.5 * sum(Bvals.^2);

%% gradient computation
Tvals = Zvals - Bvals;

G = cell(N,1);
for n = 1:N
    G{n} = zeros(size(A{n}));
end

for r = 1:R
    if (memflag)
        bigAr = cell(N,1);
        for n = 1:N
            bigAr{n} = A{n}(Wsubs(:,n),r);
        end
        for SkipN = 1:N
            newvals = Tvals;
            for n = [1:SkipN-1,SkipN+1:N]
                newvals = newvals .* bigAr{n};
            end
            G{SkipN}(:,r) = accumarray(Wsubs(:,SkipN),newvals,[sz(SkipN) 1]);
        end
    else
        for SkipN = 1:N
            newvals = Tvals;
            for n = [1:SkipN-1,SkipN+1:N]
                bigArn = A{n}(Wsubs(:,n),r);
                newvals = newvals .* bigArn;
            end
            G{SkipN}(:,r) = accumarray(Wsubs(:,SkipN),newvals,[sz(SkipN) 1]);
        end
    end

end

for n = 1:N
    G{n} = -G{n};
end
