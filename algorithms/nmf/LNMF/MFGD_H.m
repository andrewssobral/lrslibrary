function [H2,Z2,sLen2,iter]=MFGD_H(V,W,H,Z,Lpos,Lneg,tol)

% MFGD_H: multiple step length line search via Newton method for NPAF with 'active set' trick
% Model: V = WH
% Written by Naiyang Guan (ny.guan@gmail.com)
% Copyright 2010-2012 by Naiyang Guan and Dacheng Tao

% <Input>:
%       V: data matrix;
%       W: basis matrix (fixed);
%       H: current factor (target);
%       Z: current reconstruction;
%       Lpos: positive part of L;
%       Lneg: negative part of L;
%       tol: tolerence,
% <Output>:
%       H2: new value;
%       Z2: new reconstruction;
%       SLen2: step length;
%       iter: iteration number.

% Calculate scaled negative gradient at H
n=size(V,2);
G=H.*((W'*(V./Z)+H*Lneg)./(sum(W)'*ones(1,n)+H*Lpos))-H;
Z1=W*G;
C=Z./Z1;

% Hessian matrix approximation
L=Lpos-Lneg;
HSN=L.*(G'*G);
[U,S,Ut]=svds(HSN,10);

% First-order parameters
B=sum((H*L).*G)';
D=sum(Z1)';

% Newton method for multiple line search
sLen=ones(n,1);
for iter=1:20,
    [SUM1,SUM2]=SumC(V,C,sLen);
    Grad=SUM1-D-B-HSN'*sLen;    % Negative gradient
    
    % Hessian inverse
    Ainv=sparse(diag(1./SUM2));
    Hinv=inv(S+Ut'*Ainv*U);
    sLen1=sLen+(Grad-U*(Hinv*(Ut'*(Grad./SUM2))))./SUM2;
    
    if norm(sLen1-sLen)<tol,
        break;
    end
    % Update step length
    sLen=sLen1;
end

% Step length to the boundary of positive orthant
tal=0.99*ones(n,1);
if min(G(:))>=0,
   sLen3 = Inf*ones(n,1);
else
    C=H./G;
    C(C>=0)=-Inf;
    sLen3=min(-C,[],1)';        % Step size to the boundary
    tal = 0.99+0.01./sLen3;     % Modified at 4-18-2011
end

% Best step length
sLen2 = min(sLen1,tal.*sLen3);  % Modified at 4-18-2011

% Updating
ssLen2 = sparse(diag(sLen2));
H2 = H+G*ssLen2;
Z2 = Z+Z1*ssLen2;

return;