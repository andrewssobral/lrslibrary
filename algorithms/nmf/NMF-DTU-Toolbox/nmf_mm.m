function [W,H]=nmf_eucl(X,K,maxiter,speak)
%
% NMF using euclidean distance update equations :
% Lee, D..D., and Seung, H.S., (2001), 'Algorithms for Non-negative Matrix
% Factorization', Adv. Neural Info. Proc. Syst. 13, 556-562.
%
% INPUT:
% X (N,M) : N (dimensionallity) x M (samples) non negative input matrix
% K       : Number of components
% maxiter : Maximum number of iterations to run
% speak   : prints iteration count and changes in connectivity matrix 
%           elements unless speak is 0
%
% OUTPUT:
% W       : N x K matrix
% H       : K x M matrix
%
% Kasper Winther Joergensen
% Informatics and Mathematical Modelling
% Technical University of Denmark
% kwj@imm.dtu.dk
% 2006/11/16

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User adjustable parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print_iter = 50; % iterations between print on screen and convergence test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test for negative values in X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if min(min(X)) < 0
    error('Input matrix elements can not be negative');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize random W and H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n,m]=size(X);
W=rand(n,K);
H=rand(K,m);

% use W*H to test for convergence
Xr_old=W*H;

for iter=1:maxiter
    % Euclidean multiplicative method
    H = H.*(W'*X)./((W'*W)*H+eps);
    W = W.*(H*X')'./(W*(H*H')+eps);

    % print to screen
    if (rem(iter,print_iter)==0) & speak,
        Xr = W*H;
        diff = sum(sum(abs(Xr_old-Xr)));
        Xr_old = Xr;
        eucl_dist = nmf_euclidean_dist(X,W*H);
        errorx = mean(mean(abs(X-W*H)))/mean(mean(X));
        disp(['Iter = ',int2str(iter),...
            ', relative error = ',num2str(errorx),...
            ', diff = ', num2str(diff),...
            ', eucl dist ' num2str(eucl_dist)])
        if errorx < 10^(-5), break, end
    end
end
