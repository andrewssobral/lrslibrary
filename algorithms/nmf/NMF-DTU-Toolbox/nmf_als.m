function [W,H]=nmf_als(X,K,Nitsmax,speak)
%
% Truncated linear solution to LS -NMF
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
% Lars Kai Hansen, IMM-DTU (c) October 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User adjustable parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print_iter = 50; % iterations between print on screen

[D,N]=size(X);
Xscale=sum(sum(X));
%INIT
W=rand(D,K);
H=rand(K,N);
Rscale=sum(sum(W*H));
sqrnorm=sqrt(Rscale/Xscale);
H=H/sqrnorm;
W=W/sqrnorm;

Xr_old = W*H;

%ITERATE
for n=1:Nitsmax,
    %    W=X*(pinv(H*H')*H)'; % old updates
    W = ((pinv(H*H')*H)*X')';
    W=(W>0).*W;
    W=W./(repmat(sum(W),D,1)+eps); % normalize columns to unit length

    H=(W*pinv(W'*W))'*X;
    H=H.*(H>0);

    % print to screen
    if (rem(n,print_iter)==0) & speak,
        Xr = W*H;
        diff = sum(sum(abs(Xr_old-Xr)));
        Xr_old = Xr;
        eucl_dist  = nmf_euclidean_dist(X,W*H);
        errorx=mean(mean(abs(X-W*H)))/mean(mean(X));
        disp(['Iter = ',int2str(n),...
            ', relative error = ',num2str(errorx),...
            ', diff = ', num2str(diff),...
            ', eucl dist ' num2str(eucl_dist)])
        if errorx < 10^(-5), break, end
    end
end
