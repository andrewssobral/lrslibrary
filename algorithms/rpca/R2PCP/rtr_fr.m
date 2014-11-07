%% retraction of d in T_X(M_r) onto M_r
%   X=U*diag(S)*V'
%   d=U*M*V'+U_p*V'+U*V_p'
%   cf. [Van12]

function [S1,U1,V1]=rtr_fr(d,S,U,V)

n1=size(U,1);
n2=size(V,1);
n3=size(S,1);

sig=0;
if size(S,2)>1, S=diag(S); sig=1; end

eps=1e-3;

% Pu=U*U'; Pv=V*V';

M=U'*d*V;
% U_p=d*V-Pu*d*V; V_p=d'*U-Pv*d'*U;
U_p=d*V-U*M; V_p=d'*U-V*M';

[Q_u,R_u]=qr(U_p,0);
[Q_v,R_v]=qr(V_p,0);
M1=[diag(S)+M, R_v'; R_u, zeros(n3)];
[U1,S1,V1]=svd(M1);
S1=diag(S1); S1=S1(1:n3);
U1=[U,Q_u]*U1(:,1:n3);
V1=[V,Q_v]*V1(:,1:n3);

S1=max(S1,eps);

if sig==1, S1=diag(S1); end
