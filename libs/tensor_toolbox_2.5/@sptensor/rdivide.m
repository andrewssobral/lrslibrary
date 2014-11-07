function C = rdivide(A,B)
%RDIVIDE Array right division for sparse tensors.
%
%   RDIVIDE(A,B) is called for the syntax 'A ./ B' when A or B is a sparse
%   tensor. A and B must have the same size, unless one is a scalar. 
%
%   See also SPTENSOR.
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


%% Observations for sparse matrix case.
% a ./ 5 -> sparse
% 5 ./ a -> dense!
% a ./ full(a) -> sparse!
% full(a) ./ a -> dense

% Divide by a scalar -> result is sparse
if isscalar(B)
    C = mrdivide(A,B);
    return;
end

% Scalar divided by a tensor -> result is dense
if isscalar(A)
    C = A ./ full(B);
    return;
end

% Tensor divided by a tensor
if ~isequal(size(A),size(B))
    error('Must be two tensors of the same size');
end

% Two sparse tensors
if isa(A,'sptensor') && isa(B,'sptensor')

    % Find where their zeros are
    if isempty(A.subs)
        Azerosubs = allsubs(A);
    else
        Azerosubs = setdiff(allsubs(A),A.subs,'rows');
    end
    if isempty(B.subs)
        Bzerosubs = allsubs(B);
    else
        Bzerosubs = setdiff(allsubs(B),B.subs,'rows');
    end

    % Both nonzero
    [newsubs,ia,ib] = intersect(A.subs,B.subs,'rows');
    newvals = A.vals(ia) ./ B.vals(ib);
    
    % A nonzero and B zero
    moresubs = intersect(A.subs,Bzerosubs,'rows');
    morevals = repmat(Inf, size(moresubs,1),1);
    newsubs = [newsubs; moresubs];
    newvals = [newvals; morevals];
    
    % Both zero
    moresubs = intersect(Azerosubs,Bzerosubs,'rows');
    morevals = repmat(NaN, size(moresubs,1),1);
    newsubs = [newsubs; moresubs];
    newvals = [newvals; morevals];

    C = sptensor(newsubs,newvals,size(A));
    return;  
end

% Some other tensor type!
switch class(B)
  case {'tensor'}
        csubs = A.subs;
        cvals = A.vals ./ B(csubs); 
        C = sptensor(csubs, cvals, size(A));
        return;       
    case {'ktensor'}    
        R = numel(B.lambda);
        N = ndims(A);
        NZ = nnz(A);
        csubs = A.subs;
        avals = A.vals;       
        bvals = zeros(NZ,1);
        for r = 1:R
            tvals = B.lambda(r) * ones(NZ,1);
            for n = 1:N
                v = B{n}(:,r);
                tvals = tvals .* v(csubs(:,n));
            end
            bvals = bvals + tvals;
        end
        cvals = avals ./ bvals;
        C = sptensor(csubs, cvals, size(A));
        return;       
end

error('Invalid arguments for RDIVIDE.');

