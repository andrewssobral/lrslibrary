function C = times(A,B)
%TIMES Array multiplication for sparse tensors.
%
%   TIMES(A,B) is called for the syntax 'A .* B' when A or B is a 
%   sparse tensor. A and B must have the same size, unless one is a scalar.
%   A scalar can be multiplied by a sparse tensor of any size. 
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
% The result of a .* 5 is sparse.
% The result of a .* 0 is sparse.
% The result of a .* full(a) is sparse.

%%
if isscalar(B)
    C = sptensor(A.subs, A.vals * B, size(A));
    return;
end

if isscalar(A)
    C = sptensor(B.subs, B.vals * A, size(B));
    return;
end

if ~isequal(size(A),size(B))
    error('Must be two tensors of the same size');
end

switch class(B)
    case {'sptensor'}
        [csubs,ia,ib] = intersect(A.subs,B.subs,'rows');
        cvals = A.vals(ia) .* B.vals(ib);
        C = sptensor(csubs, cvals, size(A));
        return;
    case {'tensor'}
        csubs = A.subs;
        cvals = A.vals .* B(csubs); 
        C = sptensor(csubs, cvals, size(A));
        return;       
    case {'ktensor'}    
        csubs = A.subs;
        cvals = zeros(size(A.vals));       
        R = numel(B.lambda);
        N = ndims(A);
        for r = 1:R
            tvals = B.lambda(r) * A.vals;
            for n = 1:N
                v = B{n}(:,r);
                tvals = tvals .* v(csubs(:,n));
            end
            cvals = cvals + tvals;
        end
        C = sptensor(csubs, cvals, size(A));
        return;       
    otherwise
        error('Invalid second argument for sptensor/times');
end
