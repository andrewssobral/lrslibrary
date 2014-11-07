function X = normalize(X,N,normtype,mode)
%NORMALIZE Normalizes the columns of the factor matrices.
%
%   NORMALIZE(X) normalizes the columns of each factor matrix, absorbing
%   the excess weight into lambda. Also ensures that lambda is positive.  
%
%   NORMALIZE(X,N) absorbs the weights into the Nth factor matrix instead
%   of lambda. (All the lambda values are 1.)
%
%   NORMALIZE(X,0) equally divides the weights across the factor matrices.
%   (All the lambda values are 1.)
%
%   NORMALIZE(X,[]) is equivalent to NORMALIZE(X). 
%
%   NORMALIZE(X,'sort') is the same as the above except it sorts the
%   components by lambda value, from greatest to least. 
%
%   NORMALIZE(X,V,1) normalizes using the vector one norm (sum(abs(x))
%   rather than the two norm (sqrt(sum(x.^2))), where V can be any of the
%   second arguments decribed above.
%
%   NORMALIZE(X,[],1,I) just normalizes the I-th factor using whatever norm
%   is specified by the 3rd argument (1 or 2).
%
%   See also KTENSOR, NCOMPONENTS, ARRANGE.
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


%%
if ~exist('N','var')
    N = -1;
end

if isempty(N)
    N = -1;
end

if isequal(N,'sort')
    N = -2;
end

if ~exist('normtype','var')
    normtype = 2;
end

if exist('mode', 'var')
    for r = 1:length(X.lambda)
        tmp = norm(X.u{mode}(:,r),normtype);
        if (tmp > 0)
            X.u{mode}(:,r) = X.u{mode}(:,r) / tmp;
        end
        X.lambda(r) = X.lambda(r) * tmp;
    end
    return;
end

%% Ensure that matrices are normalized
for r = 1:length(X.lambda)
    for n = 1:ndims(X)
        tmp = norm(X.u{n}(:,r),normtype);
        if (tmp > 0)            
            X.u{n}(:,r) = X.u{n}(:,r) / tmp;
        end
        X.lambda(r) = X.lambda(r) * tmp;        
    end
end

%% Check that all the lambda values are positive
idx = find(X.lambda < 0);
X.u{1}(:,idx) = -1 * X.u{1}(:,idx);
X.lambda(idx) = -1 * X.lambda(idx);

%% Absorb the weight into one factor, if requested
if (N == 0)
    D = diag(nthroot(X.lambda,ndims(X)));
    X.u = cellfun(@(x) x*D, X.u, 'UniformOutput', false);
    X.lambda = ones(size(X.lambda));
elseif (N > 0)
    X.u{N} = X.u{N} * diag(X.lambda);
    X.lambda = ones(size(X.lambda));
elseif (N == -2)
    [~,p] = sort(X.lambda,'descend');
    X = arrange(X,p);
end



