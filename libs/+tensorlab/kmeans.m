function [L,C] = kmeans(X,k,dist)
%KMEANS Cluster multivariate data using the k-means++ algorithm.
%   [L,C] = kmeans(X,k) produces a 1-by-size(X,2) vector L with one class
%   label per column in X and a size(X,1)-by-k matrix C containing the
%   centers corresponding to each class.
%
%   kmeans(X,k,'norm') and kmeans(X,k,'angle') use the Euclidian distance
%   and angle, respectively, as metrics for the intra-point distance.
%
%   See also gap.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] J. B. MacQueen, "Some Methods for Classification and Analysis of 
%       MultiVariate Observations," in Proc. of the fifth Berkeley
%       Symposium on Mathematical Statistics and Probability, L. M. L. Cam
%       and J. Neyman, eds., Vol. 1, UC Press, 1967, pp. 281-297.
%   [2] D. Arthur and S. Vassilvitskii, "k-means++: The Advantages of
%       Careful Seeding," Technical Report 2006-13, Stanford InfoLab, 2006.

if nargin < 3, dist = 'norm'; end
if strcmpi(dist,'angle'), X = bsxfun(@rdivide,X,sqrt(dot(X,X,1))); end
L = [];
L1 = 0;

while length(unique(L)) ~= k
    
    % The k-means++ initialization.
    C = X(:,randi(size(X,2),1));
    L = ones(1,size(X,2));
    for i = 2:k
        D = X-C(:,L);
        D = sqrt(dot(D,D,1));
        C(:,i) = X(:,find(rand < cumsum(D)/sum(D),1));
        switch dist
            case 'norm'
                [~,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C,1).'),[],1);
            case 'angle'
                [~,L] = max(real(C'*X),[],1);
        end
    end
    
    % The k-means algorithm.
    while any(L ~= L1)
        L1 = L;
        for i = 1:k, l = L==i; C(:,i) = sum(X(:,l),2)/sum(l); end
        switch dist
            case 'norm'
                [~,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C,1).'),[],1);
            case 'angle'
                C = bsxfun(@rdivide,C,sqrt(dot(C,C,1)));
                [~,L] = max(real(C'*X),[],1);
        end
    end
    
end
