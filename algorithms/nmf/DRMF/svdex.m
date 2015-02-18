function [U s V] = svdex(X, k, options)
% [U s V] = svdex(X, k, init)
% a comprehensive svd wrapper that can determine which method 
% to use and accepts starting point
% X: the data matrix
% k: the rank of the decomposition
% options: a struct containing additional input
%   svd_solver: the svd method. can be:
%       'svd': matlab's full svd
%       'svds': matlab's sparse svd
%       'lansvd': the propack solver
%   init: the start point of the svd solver. only effective when
%         'svd_solver' is 'lansvd'.
% note that the singular values are returned as the vector s.
%
% Examples:
% [U s V] = svdex(X);
% performs the full svd of X so that X = bsxfun(@times, U, s')*V';
%
% [U s V] = svdex(X, k);
% performs the rank-k svd of X so that X ~ bsxfun(@times, U, s')*V' using the largest singular values;
% the method is selected to (approximately) maxmize speed based on X and k.
%
% % [U s V] = svdex(X, k, struct('svd_solver', 'lansvd', 'init', randn(size(X,1),1)));
% performs the rank-k svd of X using the propack solver starting from a random point.
%
% author: Liang Xiong (lxiong@cs.cmu.edu)

[m, n] = size(X);

if nargin < 3; options = []; end
[init svd_solver] = GetOptions(options, 'init', [], 'svd_solver', []);

if nargin == 1 % by default, perform full svd
    k = min(m, n);
end

if ~isempty(svd_solver)
    switch lower(svd_solver)
      case 'svd'
        [U S V] = svd(X, 'econ');
        if k < min(m,n)
            U = U(:, 1:k);
            V = V(:, 1:k);
            S = S(1:k, 1:k);
        end
      case 'svds'
        [U S V flag] = svds(X, k);
        if flag > 0 % fall back if svds failed
            warning('svds failed. falling back to svd');
            [U S V] = svd(X, 'econ');
            if k < min(m,n)
                U = U(:, 1:k);
                V = V(:, 1:k);
                S = S(1:k, 1:k);
            end
        end
      case 'lansvd'
        try
            if ~isempty(init)
                [U S V] = lansvd(X, k, 'L', struct('p0', init));
            else
                [U S V] = lansvd(X, k, 'L');
            end    
        catch ex % fall back if propack failed
            warning('lansvd failed. falling back to svds');
            [U S V] = svds(X, k);
        end
      otherwise
        error('unknown solver')
    end
    s = diag(S);
else % select solver automatically
    if issparse(X) || ChooseSVD(min(m,n),k) > 0 
        [U s V] = svdex(X, k, struct('init',init,'svd_solver','lansvd'));
    else
        [U s V] = svdex(X, k, struct('svd_solver','svd'));
    end
end

if nargout <= 1 % only need the singular values
    U = s;
end

function y = ChooseSVD(n, d)
% y = ChooseSVD( n, d)
% returns 1 if sparse svd should be used.
% This function is copied from Yi Ma's RPCA package.

if n <= 100 
    if d / n <= 0.02
        y = 1;
    else
        y = 0;
    end
elseif n <= 200
    if d / n <= 0.06
        y = 1;
    else
        y = 0;
    end
elseif n <= 300
    if d / n <= 0.26
        y = 1;
    else
        y = 0;
    end
elseif n <= 400
    if d / n <= 0.28
        y = 1;
    else
        y = 0;
    end
elseif n <= 500
    if d / n <= 0.34
        y = 1;
    else
        y = 0;
    end
else
    if d / n <= 0.38
        y = 1;
    else
        y = 0;
    end
end
