function [L,S,RMSE,error]=SSGoDec(X,rank,card,power)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Semi-Soft GoDec Algotithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS:
%X: nxp data matrix with n samples and p features
%rank: rank(L)<=rank
%card: card(S)<=card
%power: >=0, power scheme modification, increasing it lead to better
%accuracy and more time cost
%OUTPUTS:
%L:Low-rank part
%S:Sparse part
%RMSE: error
%error: ||X-L-S||/||X||
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%REFERENCE:
%Tianyi Zhou and Dacheng Tao, "GoDec: Randomized Lo-rank & Sparse Matrix
%Decomposition in Noisy Case", ICML 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tianyi Zhou, 2011, All rights reserved.

%iteration parameters
iter_max=1e+2;
error_bound=1e-3;
iter=1;
RMSE=[];

%matrix size
[m,n]=size(X);
if m<n X=X'; end

%initialization of L and S
L=X;
S=sparse(zeros(size(X)));

tic;
while true

    %Update of L
    Y2=randn(n,rank);
    for i=1:power+1
        Y1=L*Y2;
        Y2=L'*Y1;
    end
    [Q,R]=qr(Y2,0);
    L_new=(L*Q)*Q';
    
    %Update of S
    T=L-L_new+S;
    L=L_new;
    S=wthresh(T,'s',card);
    %[~,idx]=sort(abs(T(:)),'descend');
    %S=zeros(size(X));S(idx(1:card))=T(idx(1:card));
    
    %Error, stopping criteria
    %T(idx(1:card))=0;
    T=T-S;
    RMSE=[RMSE norm(T(:))];
    if RMSE(end)<error_bound || iter>iter_max
        break;
    else
        L=L+T;
    end
    iter=iter+1;
    
end
toc;

LS=L+S;
error=norm(LS(:)-X(:))/norm(X(:));
if m<n 
    LS=LS';
    L=L';
    S=S';
end