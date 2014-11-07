function [H2,Z2,sLen2,iter]=FGD(V,W,H,Z,tol)

% Fast gradient descent (FGD) with Newton method for NNLS
% Copyright 2010-2012 by Naiyang Guan and Dacheng Tao
% Arguments:
%     sLen: step length

% Calculate scaled negative gradient at W
n=size(V,2);
G=H.*(W'*(V./Z))./(sum(W)'*ones(1,n))-H;
Z1=W*G;
d=sum(sum(Z1));
C=Z./Z1;

% Newton method
sLen=1;
for iter=1:20
    [sum1,sum2]=SumC1(V,C,sLen);
    sLen1=sLen-(d-sum1)/sum2;
    if abs(sLen1-sLen)<tol
        break;
    else
        sLen=sLen1;
    end
end

% Step length to boundary of positive orthant
if min(G(:))>=0
    sLen3=Inf;
else
    C=H./G;
    sLen3=min(-C(C<0));
end

% Best step length
sLen2=max(min(sLen1,0.99*sLen3),1);
H2=H+sLen2*G;
Z2=Z+sLen2*Z1;

return;