function [ Z, trnorm, U ] = matrix_shrinkage( X, lambda )
% X is a matrix
% Shrinkage operation on a matrix X = U*diag(s)*V'
% Z = U*diag(max(0,s-lambda))*V'
% opt = 'matrix' | 'tensor'

% uncomment the code to use PROPACK to compute partial SVD

global sv;
global tmode;
global use_propack;

n = min(size(X,1),size(X,2));
i = tmode;
sv_local = sv(i);

if choosvd( n, sv_local) && use_propack
    opt.delta = eps;
    opt.eta = eps;
    [ U, S, V ] = lansvd( X, sv_local, 'L', opt);
else
    [ U, S, V ] = svd( X, 'econ' );     %fprintf( 'sv(%d) = %d \n', i, sv_local );
end
s = diag(S);
s = s - lambda;     sMax = max(s);
posI = s > 0;   %sum(posI)
Z = scale_matrix( U(:,posI), s(posI), 1 ) * V(:,posI)';
trnorm = sum(s(posI));

% adjust sv
% svp should be no. of sv's larger than a threshold!!
svp = sum(posI);
if svp < sv_local
    sv_local = min(svp + 1, n);
else
    sv_local = min(svp + round(0.05*n), n);
end
sv(i) = sv_local;
end