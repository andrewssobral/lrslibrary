function K = fixsigns(K,K0)
%FIXSIGNS Fix sign ambiguity of a ktensor.
%
%   K = FIXSIGNS(K) makes it so that the largest magnitude entries for
%   each vector in each factor of K are positive, provided that the
%   sign on *pairs* of vectors in a rank-1 component can be flipped.
%
%   K = FIXSIGNS(K,K0) returns a version of K where some of the signs of
%   the columns of the factor matrices have been flipped to better align
%   with K0. 
%
%   See also KTENSOR and KTENSOR/ARRANGE.
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


if nargin == 1
    K = fixsigns_oneargin(K);
else 
    K = fixsigns_twoargin(K, K0);
end


%%
function K = fixsigns_oneargin(K)
R = length(K.lambda);
for r = 1 : R
    
    for n = 1:ndims(K)
        [val(n),idx(n)] = max(abs(K.u{n}(:,r)));    
        sgn(n) = sign(K.u{n}(idx(n),r));
    end

    negidx = find(sgn == -1);
    nflip = 2 * floor(numel(negidx)/2);

    for i = 1:nflip
        n = negidx(i);
        K.u{n}(:,r) =  -K.u{n}(:,r);
    end

end

%%
function [A, best_sign] = fixsigns_twoargin(A,B)
%   T. Kolda, November 2010.

if ~isa(A,'ktensor')
    A = ktensor(A);
end
if ~isa(B,'ktensor')
    B = ktensor(B);
end
A = normalize(A);
B = normalize(B);

N = ndims(A);
RA = ncomponents(A);
RB = ncomponents(B);

%% Try to fix the signs for each component
best_sign = ones(N,RA);
for r = 1:RB
    
    % Compute the inner products. They should mostly be O(1) if there is a
    % good match because the factors have prevsiouly been normalized. If
    % the signs are correct, then the score should be +1. Otherwise we need
    % to flip the sign and the score should be -1.
    sgn_score = zeros(N,1);
    for n = 1:N
        sgn_score(n) = A{n}(:,r)' * B{n}(:,r);
    end
    
    % Sort the sign scores.
    [sort_sgn_score, sort_idx] = sort(sgn_score,'ascend');
    
    % Determine the number of scores that should be flipped.
    breakpt = find(sort_sgn_score < 0, 1, 'last');
    
    % If nothing needs to be flipped, then move on the the next component.
    if isempty(breakpt)
        continue;
    end
    
    % Need to flip signs in pairs. If we don't have an even number of
    % negative sign scores, then we need to decide to do one fewer or one
    % more.
    if (mod(breakpt,2) == 0)
        endpt = breakpt;
    else
        fprintf('Trouble fixing signs for mode %d\n', r);
        if (breakpt < RB) && (-sort_sgn_score(breakpt) > sort_sgn_score(breakpt+1))
            endpt = breakpt + 1;
        else
            endpt = breakpt - 1;
        end        
    end    
    
    % Flip the signs
    for i = 1:endpt
        A{sort_idx(i)}(:,r) = -1 * A{sort_idx(i)}(:,r);
        best_sign(sort_idx(i),r) = -1;
    end
       
end



