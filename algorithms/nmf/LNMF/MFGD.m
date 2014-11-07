function [H2,Z2,sLen2,iter]=MFGD(V,W,H,Z,tol)

% MFGD: multiple step length line search via Newton method for NNLS with 'active set' trick
% Model: V = WH
% Written by Naiyang Guan (ny.guan@gmail.com)
% Copyright 2010-2012 by Naiyang Guan and Dacheng Tao

% <Input>:
%       V: data matrix;
%       W: basis matrix (fixed);
%       H: current factor (target);
%       Z: current reconstruction;
%       tol: tolerence,
% <Output>:
%       H2: new value;
%       Z2: new reconstruction;
%       SLen2: step length;
%       iter: iteration number.

% Calculate scaled negative gradient at H
n=size(V,2);
G=H.*((W'*(V./Z))./(sum(W)'*ones(1,n)))-H;
Z1=W*G;
d=sum(Z1)';
C=Z./Z1;

% Newton method for multiple line search
sLen=ones(n,1);
index=ones(1,n);    % Free set index
for iter=1:20,
    iindex=(index==1);
    [sum1,sum2]=SumC(V(:,iindex),C(:,iindex),sLen(iindex));
    sLen1=sLen;
    sLen1(iindex)=sLen(iindex)+(sum1-d(iindex))./sum2;
    index(abs(sLen1-sLen)<tol)=0;
    if sum(index)==0
        break;
    end
    sLen=sLen1;
end

% Step length to boundary of positive orthant
tal=0.99*ones(n,1);
if min(G(:))>=0
    sLen3=Inf*ones(n,1);
else
    C=H./G;
    C(C>=0)=-Inf;
    sLen3=min(-C,[],1)';
    tal=0.99+0.01./sLen3;
end

% Best step length
sLen2=min(sLen1,tal.*sLen3);
ssLen2=sparse(diag(sLen2));
H2=H+G*ssLen2;
Z2=Z+Z1*ssLen2;

return;