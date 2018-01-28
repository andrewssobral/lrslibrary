function P_hat = simpleEVD(X, r)

%%%Function to implement the Simple-EVD algorithm and returns a basis
%%%matrix for the new subspace -- edit these comments before Gitting

%%%                         Inputs                          %%%
%%%         X - Emperical Covariance data matrix (m X m)    %%%
%%%         r - Target rank of output                       %%%

%%%                         Outputs                         %%%
%%%         P_hat - Basis matrix for output (m X r)         %%%

[P_hat, ~] = svds(X, r);
end