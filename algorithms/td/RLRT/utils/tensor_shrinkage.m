function [ T, trnorm, U ] = tensor_shrinkage( X, tau, mode )
% X and T are tensors
global tmode;
tmode = mode;

Xmat = tenmat(X,mode);
[Tmat, trnorm, U] = matrix_shrinkage( double(Xmat), tau );
% [ U, S, V ] = svd( double(Xmat), 'econ' );
% s = diag(S);
% s = s - tau;
% s( s < 0 ) = 0;
% Tmat = scale_matrix( U, s, 1 ) * V';
T = tensor( tenmat( Tmat, Xmat.rdims, Xmat.cdims, Xmat.tsize ) );

end