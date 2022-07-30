function Y = optimal_mc_ls_Y( W, b, Omega, Zs, N, grad_sum, mu, lambda )

R = tenzeros( size(W) );
R(Omega) = lambda*b;
R = R - (grad_sum - Zs{N+1}) + (N/mu)*W;
Y = R * mu/N;
Y(Omega) = R(Omega) / (lambda+N/mu);

end