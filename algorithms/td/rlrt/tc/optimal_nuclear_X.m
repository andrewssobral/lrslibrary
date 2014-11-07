function [ X, U, gammax, V ] = optimal_nuclear_X( W, grad_sum, Zs, i, N, mu, sig )

temptenmat = tenmat( N*W - mu*(grad_sum-Zs{i}), i );
[ U, S, V ] = svd( double( temptenmat ), 'econ' );
eta = diag(S);
gammax = zeros( size(eta) );
ind1 = eta <= sig*N + mu;   %sum(ind1==1)
gammax(ind1) = sig*eta(ind1) / (sig*N+mu);
gammax(~ind1) = (eta(~ind1)-mu) / N;
X_i = scale_matrix( U, gammax, 1 ) * V';
X = tensor( tenmat(X_i, temptenmat.rdims, temptenmat.cdims, temptenmat.tsize) );
        
end