
function [U1,S]=BL_EVD(A,U0,k)

% Jan 2010
% This matlab code implements the block lanczos method with k steps
% for partial EVD.
%
% A - m * m symmetric matrix (required input)
% U0 - Initial subspace (required input)
% k - the number of blocks
%
% U1 - Output subspace, eigenvectors
% S - eigenvalues
%
%
% Reference: Zhouchen Lin and Siming Wei, A Block Lanczos with Warm Start 
% Technique for Accelerating Nuclear Norm Minimization Algorithms, 
% arxiv: 1012.0365.
%
% Bug report: zclin2000@hotmail.com
%


if nargin < 2
    error('Too few arguments') ;
elseif size(A,1)~=size(A, 2) || norm(A-A','fro')>1e-8*norm(A,'fro')
    error('A should be symmetric!');
end

if nargin < 3
    k = 4;
elseif k == -1
    k = 4;
end

[m,r]=size(U0);
% tridiagolization;    
X(:,:,1)=U0;
AX=A*U0;
M(:,:,1)=U0'*AX;
Q=X(:,:,1);
Bdiag=[];
Mdiag=M(:,:,1);
l=1;
stop=0;

while stop~=1
    if l==1
    R(:,:,l)=AX-U0*M(:,:,1);
    else
    R(:,:,l)=AX-X(:,:,l)*M(:,:,l)-X(:,:,l-1)*B(:,:,l-1)';
    end
    [X(:,:,l+1),B(:,:,l)]=qr(R(:,:,l),0); 
    AX=A*X(:,:,l+1);
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

% Now T=Q'CQ is block tridiagonal;
B=[zeros(r,l*r);Bdiag,zeros(l*r-r,r)];
T=Mdiag+B+B';
T=(T+T')/2;

% Compute EVD of T: T=Z*S*Z';
[Z,S0]=eig(T);
diagS=diag(S0);
S=diag(diagS(l*r:-1:l*r-r+1));
U1=Q*Z(:,l*r:-1:l*r-r+1);
