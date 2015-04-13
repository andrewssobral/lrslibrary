function x = prox_matrix(v, lambda, prox_f)
% PROX_MATRIX    The proximal operator of a matrix function.
%
%   Suppose F is a orthogonally invariant matrix function such that
%   F(X) = f(s(X)), where s is the singular value map and f is some
%   absolutely symmetric function. Then
%
%     X = prox_matrix(V,lambda,prox_f)
%
%   evaluates the proximal operator of F via the proximal operator
%   of f. Here, it must be possible to evaluate prox_f as prox_f(v,lambda).
%
%   For example,
%
%     prox_matrix(V,lambda,prox_l1)
%
%   evaluates the proximal operator of the nuclear norm at V
%   (i.e., the singular value thresholding operator).

    %[U,S,V] = svd(v,'econ');
    [U,S,V] = svdecon(v); % fastest
    x = U*diag(prox_f(diag(S), lambda))*V';
end
