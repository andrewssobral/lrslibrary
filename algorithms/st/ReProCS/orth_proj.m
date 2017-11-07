function[phi] = orth_proj(P)
%%% Funciton to compute the projection of a vector orthogonal to the
%%% subspace spanned by a matrix 

%%%         Inputs
%%%     P - basis matrix
%%%     x - input vector

%%%         Outputs
%%%     phi- orthogonal projection matrix
%%%     y - output vector

[nr, ~] = size(P);

phi = eye(nr) - P * P';
%y = x - P * (P' * x);

end