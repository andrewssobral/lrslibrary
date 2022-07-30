function [ T, trnorm, U ] = tensor_hard_thresh( X, k, mode )
% X is a tensor
% Computes the k largest singular values and assoc. singular vectors for
% mode (mode)

% global sv;
% global tmode;
global use_propack;

Xmat = tenmat(X,mode);
dXmat = double( Xmat );

n = min(size(dXmat,1),size(dXmat,2));
% i = tmode;
% sv_local = sv(i);

if choosvd( n, k) && use_propack
    opt.delta = eps;
    opt.eta = eps;
    [ U, S, V ] = lansvd( dXmat, k, 'L', opt);
else
    [ U, S, V ] = svds( dXmat, k );     %fprintf( 'sv(%d) = %d \n', i, sv_local );
end
s = diag(S);
trnorm = sum(s);
Tmat = scale_matrix( U, s, 1 ) * V';

T = tensor( tenmat( Tmat, Xmat.rdims, Xmat.cdims, Xmat.tsize ) );
end