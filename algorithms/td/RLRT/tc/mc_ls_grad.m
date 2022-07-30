function grad =  mc_ls_grad( X, b, Omega, lambda, opt )
% grad = lambda*( A'A(X)+A'*b)
% X = a tensor or matrix
% b = a vector of known entries
% Omega = set of matching indices: 
%   for 'tensor': a matrix of indices (each row is an index)
%   for 'matrix': a vector of linear indices
% opt = 'tensor', 'matrix'


if exist( 'opt', 'var' ) && strcmp( opt, 'tensor' )
    grad = tenzeros( size( X ) );
    
else
    grad = zeros( size( X ) );    
end

grad( Omega ) = lambda *( X( Omega ) - b ); 


end