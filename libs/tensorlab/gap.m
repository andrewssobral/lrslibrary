function [L,C,gap,s] = gap(X,options)
%GAP Optimal clustering based on the gap statistic.
%   gap(X) plots the gap statistic [1] of different clusterings of the
%   dataset X, in which the rows are variables and the columns are
%   observations, and marks the optimal number of clusters with a square.
%
%   [L,C] = gap(X) does not plot anything and produces a 1-by-size(X,2)
%   vector L with one class label per column in X and a size(X,1)-by-k
%   matrix C containing the centers corresponding to each class, where k is
%   the optimal number of clusters as determined by the gap statistic. The
%   algorithm stops as soon as k is determined.
%
%   [L,C,gap,s] = gap(X) additionally returns the gap statistic gap and
%   standard deviations s used to create the gap statistic plot above.
%
%   gap(X,options) may be used to set the following options:
%
%      options.B = 8            - The number of reference datasets to
%                                 compute the dispersion of with kmeans.
%      options.Dist = 'norm'    - The distance metric to use in the kmeans
%                                 algorithm, see help kmeans.
%      options.MaxK = size(X,2) - The maximal number of clusters.
%      options.MinK = 1         - The minimal number of clusters.
%      options.Tries = 8        - For each number of clusters, kmeans is
%                                 executed options.Tries times on X and
%                                 the result with the smallest dispersion
%                                 is kept.
%
%   See also kmeans.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] R. Tibshirani, G. Walther, T. Hastie, "Estimating the number of
%       clusters in a data set via the gap statistic," J. R. Statist. Soc.
%       B, Vol. 63, No. 2, 2001, pp. 411-423.
%   [2] D. E. Knuth, "The Art of Computer Programming, vol. 2:
%       Seminumerical Algorithms," Addison-Wesley, 1969.

% Check the options structure.
if nargin < 2, options = struct; end
if ~isfield(options,'B'), options.B = 5; end
if ~isfield(options,'Dist'), options.Dist = 'norm'; end
if ~isfield(options,'MaxK'), options.MaxK = size(X,2); end
if ~isfield(options,'MinK'), options.MinK = 1; end
if ~isfield(options,'Tries'), options.Tries = 5; end
if strcmpi(options.Dist,'angle')
    X = bsxfun(@rdivide,X,sqrt(dot(X,X)));
end

% Generate reference distribution.
switch options.Dist
    case 'norm'
        % Uniform points in the principal component box.
        C = mean(X,2);
        [U,S,V] = svd(bsxfun(@minus,X,C),'econ');
        mn = diag(S).*min(V).';
        mx = diag(S).*max(V).';
        Z = arrayfun(@(i)bsxfun(@plus,U*bsxfun(@plus,bsxfun(@times, ...
                rand(size(X)),mx-mn),mn),C), ...
                1:options.B,'UniformOutput',false);
    case 'angle'
        % Uniform points on a sphere [2] for lack of spherical PCA.
        Z = cell(1,options.B);
        for i = 1:options.B
            Z{i} = randn(size(X));
            Z{i} = bsxfun(@rdivide,Z{i},sqrt(dot(Z{i},Z{i})));
        end
end

% Compute gap statistic.
gap = zeros(1,options.MaxK);
s = zeros(1,options.MaxK);
kopt = nan;
for k = options.MinK:options.MaxK
    if k == 1
        Lk = ones(1,size(X,2)); Ck = mean(X,2);
        WX = dispersion(X,Lk);
        WZ = cellfun(@(z)dispersion(z,Lk),Z);
    else
        if k > options.MinK, L = Lk; C = Ck; end
        [Lk,Ck] = kmeans(X,k,options.Dist);
        WX = dispersion(X,Lk);
        for n = 2:options.Tries
            [Ln,Cn] = kmeans(X,k,options.Dist);
            WXn = dispersion(X,Ln);
            if WXn < WX, Lk = Ln; Ck = Cn; WX = WXn; end
        end
        WZ = cellfun(@(z)dispersion(z,kmeans(z,k,options.Dist)),Z);
    end
    gap(k) = mean(log(WZ))-log(WX);
    s(k) = std(log(WZ))*sqrt(1+1/options.B);
    if k > options.MinK && gap(k-1) >= gap(k)-s(k)
        if isnan(kopt), kopt = k-1; end
        if nargout == 1 || nargout == 2
            gap = gap(1:k);
            s = s(1:k);
            break;
        end
    end
end
if isnan(kopt)
    L = 1:size(X,2);
    C = X;
end
gap = gap(options.MinK:end);
s = s(options.MinK:end);

% Plot gap statistic if no output is captured.
if nargout == 0
    plot(kopt,gap(kopt-options.MinK+1),'s'); hold on;
    errorbar(options.MinK:k,gap,s); hold off;
    xlabel('k');
    ylabel('gap');
    legend(['Optimal k = ' int2str(kopt)],'Gap statistic','Location','NW');
end

function W = dispersion(X,L)
    W = 0;
    for l = 1:max(L)
        Xl = X(:,L == l);
        switch options.Dist
            case 'norm'
                D = bsxfun(@minus,permute(Xl,[2 3 1]),permute(Xl,[3 2 1]));
                W = W+(D(:)'*D(:))/(2*size(Xl,2));
            case 'angle'
                D = real(acos(real(Xl'*Xl)));
                W = W+sum(D(:))/(2*size(Xl,2));
        end
    end
end

end
