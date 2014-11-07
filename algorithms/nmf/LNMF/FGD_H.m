function [H2,Z2,sLen2,iter]=FGD_H(V,W,H,Z,Lpos,Lneg,tol)

% Fast gradient descent (FGD) with Newton method for NPAF
% Copyright@Naiyang Guan and Dacheng Tao
% Arguments:
%     sLen: step length

% Calculate scaled negative gradient at H
n=size(H,2);
G=H.*(W'*(V./Z)+H*Lneg)./(sum(W)'*ones(1,n)+H*Lpos)-H;
L=Lpos-Lneg;
Z1=W*G;
a=sum(sum(L.*(G'*G)));
d=sum(sum(L.*(G'*H)))+sum(sum(Z1));
C=Z./Z1;

% Newton method
sLen=1;
for iter=1:20,
    [sum1,sum2]=SumC1(V,C,sLen);
    sLen1=sLen-(a*sLen-sum1+d)/(a+sum2);
    if abs(sLen1-sLen)<tol,
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