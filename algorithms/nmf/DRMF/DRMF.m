function [ L, S ] = DRMF(M, K, E, options)
% [ L, S ] = DRMF(M, K, E, options)
% direct robust matrix factorization
% M: input matrix
% K: max rank
% E: max percentage of outliers
% 
% options: 
%     init: the initial L. default the raw SVD.
%     type: 'E' for element outliers, 'R' for row outliers, 'C' for column outliers. default 'E'.
%     max_iter: max iteration. default 100.
%     epsilon: relative change of objective to stop the iteration. default 1e-4.
%     verbose: show info or not
%
% output:
% L: Low rank result
% S: sparse outliers
%
% author: Liang Xiong (lxiong@cs.cmu.edu)

if nargin < 4;    options = [];    end
[init, type, max_iter, epsilon, verbose] = GetOptions(options, ...
    'init', [], 'type', 'E', 'max_iter', 100, 'epsilon', 1e-4, 'verbose', true);
[m, n] = size(M);

if isempty(init)
    % default start point
    [U, s, V] = svdex(M, K);
    L = bsxfun(@times, U, s')*V';
else
    L = init;
    [U, dummy] = svdex(L, 1, struct('svd_solver','svds'));
end

if verbose; 
    fprintf('DRMF-%s for %dx%d matrix with rank %d and %0.1f%% outliers\n', ...
            type, m, n, K, E*100);
end

tic;
objs = nan(max_iter, 1);
S = zeros(m, n);
for iter = 1:max_iter
    old = {L, S};
    
    % update S
    A = M - L;
    switch lower(type)
        case 'e'
            idx = abs(A) >= qt(abs(A), 1 - E);
            [ii, jj] = find(idx);
            S = sparse(ii, jj, A(idx), m, n);
        case 'c'
            sA = sos(A, 1);
            ii = sA >= qt(sA, 1 - E);
            S(:) = 0; S(:, ii) = A(:, ii);
        case 'r'
            sA = sos(A, 2);
            ii = sA >= qt(sA, 1 - E);
            S(:) = 0; S(ii, :) = A(ii, :);
        otherwise
            error('unknown outlier type')
    end
    
    % update L
    B = M - S;
    % you can use approximate svd here to get faster
    [U, ss, V] = svdex(B, K, struct('init',U(:,1),'svd_solver','lansvd'));
    L = bsxfun(@times, U, ss')*V';
    
    objs(iter, :) = RMSE(B - L);
    if verbose
        fprintf('-Iter=%d, Obj=%g, Time=%0.3f\n', iter, objs(iter), toc);
    end

    % check convergence
    if iter > 1 && ((objs(iter-1) - objs(iter))/objs(iter-1) < epsilon || objs(iter) < 1e-7)
        [L, S] = deal(old{:});
        break
    end
end

function q = qt(x, p)
% compute the quantile

assert(p > 0 & p < 1, 'wrong quantile');

x = sort(x(:));
q = x(max(1, round(p*length(x))));
