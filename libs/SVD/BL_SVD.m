function [U,S,V]=BL_SVD(A,u0,v0,k)

% Block lanczos method to compute SVD of matrix A with k steps; 
% Suppose A is m*n, where m>n; A0=U0*S0*V0', where U0 is m*r, V0 is n*r;
% This is equivalent to compute block lanczos on C=[0,A;A',0];
% The initial lanczos block should set to be uv0=[U0;V0]/sqrt(2);
% where r is the rank of A0 and c is the slack positive interger variable;
%
% A - m * n matrix (required input)
% [u0; v0]/sqrt(2) - Initial subspace (required input)
% k - the number of blocks
%
% U S V - partial singular value decomposition
%
% Reference: Zhouchen Lin and Siming Wei, A Block Lanczos with Warm Start 
% Technique for Accelerating Nuclear Norm Minimization Algorithms, 
% arxiv: 1012.0365.
%
% Bug report: zclin2000@hotmail.com
%

if nargin < 3
    error('Too few arguments') ;
end

if nargin < 4
    k = 1;
elseif k == -1
    k = 1;
end


[m,n]=size(A);
mark=0;
if m<n % Let m>=n;
    mark=1;
    A=A';
    t=m;
    m=n;
    n=t;
    uv0=[v0; u0]/sqrt(2);
else
    uv0=[u0; v0]/sqrt(2);
end
r=size(uv0,2);
% Initialization;

% tridiagolization;    
X(:,:,1)=uv0;
% AX=multiC(A,uv0,m,n);
AX = [A*uv0(m+1:m+n,:);A'*uv0(1:m,:)];
M(:,:,1)=uv0'*AX;
Q=X(:,:,1);
Bdiag=[];
Mdiag=M(:,:,1);
l=1;
stop=0;
tol=1.e-3;

while stop~=1
    if l==1
    R(:,:,l)=AX-uv0*M(:,:,1);
    else
    R(:,:,l)=AX-X(:,:,l)*M(:,:,l)-X(:,:,l-1)*B(:,:,l-1)';
    end
    R_max=max(max(abs(R(:,:,l))));
    
    if R_max<= tol
        stop=1;
    else
    [X(:,:,l+1),B(:,:,l)]=qr(R(:,:,l),0); 
%     AX=multiC(A,X(:,:,l+1),m,n);
%     Y(:,:) = X(:,:,l+1);
    AX = [A*X(m+1:m+n,:,l+1);A'*X(1:m,:,l+1)];
    M(:,:,l+1)=X(:,:,l+1)'*AX;
    Q=[Q,X(:,:,l+1)];
    
    % Record the tridiagonal elments;
    Bdiag=blkdiag(Bdiag,B(:,:,l));
    Mdiag=blkdiag(Mdiag,M(:,:,l+1));
    l=l+1;
         if l>=k
             stop=1;
         end
    end
end

% Now T=Q'CQ is block tridiagonal;
B0=[zeros(r,l*r);Bdiag,zeros(l*r-r,r)];
T=Mdiag+B0+B0';
T=(T+T')/2;

% Compute EVD of T: T=Z*S*Z';
[Z,S]=eig(T);

% Therefore (QZ)'*C*(QZ)=S
add=0;
QZ=Q*Z(:,l*r:-1:l*r-r+1-add);
U_temp=QZ(1:m,:)*sqrt(2);
V_temp=QZ(m+1:m+n,:)*sqrt(2);
dS=diag(S);
S=diag(dS(l*r:-1:l*r+1-r-add));
if mark==1
    U=V_temp;
    V=U_temp;
else
    U=U_temp;
    V=V_temp;
end