function [X S Y dist] = OptSpace(M_E,r,niter,tol)
% An algorithm for Matrix Reconstruction from a partially revealed set. 
% See "Matrix Completion from a Few Entries"(http://arxiv.org/pdf/0901.3150) for details
% Usage :
% [X S Y dist] = OptSpace(A,r,niter,tol);
% [X S Y dist] = OptSpace(A);
% 
% INPUT :
% A     :  The partially revealed matrix.
%          Sparse matrix with zeroes at the unrevealed indices.
%
% r     :  The rank to be used for reconstruction. Use [] to guess the rank.
% niter :  The max. no. of iterations. Use [] to use default (50).
% tol   :  Stop iterations if norm( (XSY' - M_E).*E , 'fro' )/sqrt(|E|) < tol, where
%        - E_{ij} = 1 if M_{ij} is revealed and zero otherwise, 
%        - |E| is the size of the revealed set.				
%        - Use [] to use the default (1e-6)
%
%
% OUTPUT :
% X      : A size(A,1)xr matrix
% S      : An rxr matrix
% Y      : A size(A,2)xr matrix
% such that M_hat = X*S*Y' 
% dist   : A vector containing norm( (XSY' - M_E).*E , 'fro' )/sqrt(|E|) at each
%          successive iteration
%
% Date : 21st April, 2009
% COPYRIGHT 2009 Raghunandan H. Keshavan, Andrea Montanari, Sewoong Oh





if(nargin==1)
	
	M_E = sparse(M_E);
	[n m] = size(M_E);
	E = spones(M_E);
	eps = nnz(E)/sqrt(m*n) ;

	tol = 1e-6;

	fprintf(1,'Rank not specified. Trying to guess ...\n');
	r = guessRank(M_E) ;
	fprintf(1,'Using Rank : %d\n',r);
	
	m0 = 10000 ;
	rho = 0;
	
	niter = 50;
elseif(nargin==4)
	
	M_E = sparse(M_E);
	[n m] = size(M_E);
	E = spones(M_E);
	eps = nnz(E)/sqrt(m*n) ;

	if( length(tol) == 0 )
		tol = 1e-6;
	end

	if( length(r) == 0 )

		fprintf(1,'Rank not specified. Trying to guess ...\n');
		r = guessRank(M_E) ;
		fprintf(1,'Using Rank : %d\n',r);
	end

	m0 = 10000 ;
	rho = 0;

	if( length(niter) == 0 )
		niter = 50 ;
	end	
else
	fprintf(1,'Improper arguments (See "help OptSpace")\n');
	fprintf(1,'Usage :\n[X S Y dist] = OptSpace(A,r,niter,tol) \n') ;
	fprintf(1,'[X S Y dist] = OptSpace(A)\n');
	return;
end	

rescal_param = sqrt( nnz(E) * r / norm(M_E,'fro')^2 ) ;
M_E = M_E * rescal_param ;

fprintf(1,'Trimming ...\n');
% Trimming

M_Et = M_E ;
d = sum(E);
d_=mean(full(d));
for col=1:m
    if ( sum(E(:,col))>2*d_ )
        list = find( E(:,col) > 0 );
        p = randperm(length(list));
        M_Et( list( p(ceil(2*d_):end) ) , col ) = 0;
    end
end

d = sum(E');
d_= mean(full(d));
for row=1:n
    if ( sum(E(row,:))>2*d_ )
        list = find( E(row,:) > 0 );
        p = randperm(length(list));
        M_Et(row,list( p(ceil(2*d_):end) ) ) = 0;
    end
end

fprintf(1,'Sparse SVD ...\n');
% Sparse SVD
[X0 S0 Y0] = svds(M_Et,r) ;

clear M_Et;

% Initial Guess
X0 = X0*sqrt(n) ; Y0 = Y0*sqrt(m) ;
S0 = S0 / eps ;


fprintf(1,'Iteration\tFit Error\n');

% Gradient Descent
X = X0;Y=Y0;
S = getoptS(X,Y,M_E,E);


dist(1) = norm( (M_E - X*S*Y').*E ,'fro')/sqrt(nnz(E) )  ;
fprintf(1,'0\t\t%e\n',dist(1) ) ;

for i = 1:niter

% Compute the Gradient 
	[W Z] = gradF_t(X,Y,S,M_E,E,m0,rho);

% Line search for the optimum jump length	
	t = getoptT(X,W,Y,Z,S,M_E,E,m0,rho) ;
	X = X + t*W;Y = Y + t*Z;S = getoptS(X,Y,M_E,E) ;
	
% Compute the distortion	
	dist(i+1) = norm( (M_E - X*S*Y').*E,'fro' )/sqrt(nnz(E));
	fprintf(1,'%d\t\t%e\n',i,dist(i+1) ) ;
	if( dist(i+1) < tol )
		break ;
	end
end

S = S /rescal_param ;

% Function to Guess the Rank of the input Matrix
function r = guessRank(M_E);
	[n m] = size(M_E);
	epsilon = nnz(M_E)/sqrt(m*n);
    S0 = svds(M_E,100) ;

    S1=S0(1:end-1)-S0(2:end);
    S1_ = S1./mean(S1(end-10:end));
    r1=0;
    lam=0.05;
    while(r1<=0)
        for idx=1:length(S1_)
            cost(idx) = lam*max(S1_(idx:end)) + idx;
        end
        [v2 i2] = min(cost);
        r1 = max(i2-1);
        lam=lam+0.05;
    end

	clear cost;
    for idx=1:length(S0)-1
        cost(idx) = (S0(idx+1)+sqrt(idx*epsilon)*S0(1)/epsilon  )/S0(idx);
    end
    [v2 i2] = min(cost);
    r2 = max(i2);

	r = max([r1 r2]);






% * * * * * * * * * * * * * * * * * * * *
% Function to compute the distortion
function out = F_t(X,Y,S,M_E,E,m0,rho)
[n r] = size(X) ;

out1 = sum( sum( ( (X*S*Y' - M_E).*E ).^2 ) )/2 ;

out2 =  rho*G(Y,m0,r) ;
out3 =  rho*G(X,m0,r) ;
out = out1+out2+out3 ;

function out = G(X,m0,r)

z = sum(X.^2,2)/(2*m0*r) ;
y = exp( (z-1).^2 ) - 1 ;
y( find(z < 1) ) = 0 ;
out = sum(y) ;
% * * * * * * * * * * * * * * * * * * * *



% Function to compute the gradient
function [W Z] = gradF_t(X,Y,S,M_E,E,m0,rho)
[n r] = size(X);
[m r] = size(Y);

XS = X*S ;
YS = Y*S' ;
XSY = XS*Y' ;

Qx = X'* ( (M_E - XSY).*E )*YS /n;
Qy = Y'* ( (M_E - XSY).*E )'*XS /m;

W = ( (XSY - M_E).*E )*YS + X*Qx + rho*Gp(X,m0,r);
Z = ( (XSY - M_E).*E )'*XS + Y*Qy + rho*Gp(Y,m0,r);

function out = Gp(X,m0,r)
z = sum(X.^2,2) /(2*m0*r) ;
z = 2*exp( (z-1).^2 ).*(z-1) ;
z( find(z<0) ) = 0;

out = X.*repmat(z,1,r) / (m0*r) ;
% * * * * * * * * * * * * * * * * * * * *


% * * * * * * * * * * * * * * * * * * * *
% Function to find Sopt given X, Y
function out = getoptS(X,Y,M_E,E)

[n r] = size(X);
C = X' * ( M_E ) * Y ; C = C(:) ;

for i = 1:r
        for j = 1:r
                ind = (j-1)*r + i ;
                temp = X' * (  (X(:,i) * Y(:,j)').*E ) * Y ;
                A(:,ind) = temp(:) ;
        end
end

S = A\C ;
out = reshape(S,r,r) ;
% * * * * * * * * * * * * * * * * * * * *


% * * * * * * * * * * * * * * * * * * * *
% Function to perform line search
function out = getoptT(X,W,Y,Z,S,M_E,E,m0,rho)
norm2WZ = norm(W,'fro')^2 + norm(Z,'fro')^2;
f(1) = F_t(X, Y,S,M_E,E,m0,rho) ;

t = -1e-1 ;
for i = 1:20
        f(i+1) = F_t(X+t*W,Y+t*Z,S,M_E,E,m0,rho) ;

        if( f(i+1) - f(1) <= .5*(t)*norm2WZ )
            out = t ;
            return;
        end
        t = t/2 ;
end
out = t ;


