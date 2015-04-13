function [hessfd, storedb] = getHessianFD(problem, x, d, storedb)
% Computes an approx. of the Hessian w/ finite differences of the gradient.
%
% function [hessfd, storedb] = getHessianFD(problem, x, d, storedb)
%
% Return a finite difference approximation of the Hessian at x along d of
% the cost function described in the problem structure. The cache database
% storedb is passed along, possibly modified and returned in the process.
% The finite difference is based on computations of the gradient. 
%
% If the gradient cannot be computed, an exception is thrown.

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 

    
    if ~canGetGradient(problem)
        up = MException('manopt:getHessianFD:nogradient', ...
            'getHessianFD requires the gradient to be computable.');
        throw(up);
    end
    
    % First, check whether the step d is not too small
    if problem.M.norm(x, d) < eps
        hessfd = problem.M.zerovec(x);
        return;
    end
    
    % Parameter: how far do we look?
	% TODO: this parameter should be tunable by the users.
    epsilon = 1e-4;
        
    % TODO: to make this approximation of the Hessian radially linear, that
    % is, to have that HessianFD(x, a*d) = a*HessianFD(x, d) for all
    % points x, tangent vectors d and real scalars a, we need to pay
    % special attention to the case of a < 0. This requires a notion of
    % "sign" for vectors d.
	% If vectors are uniquely represented by a matrix (which is the case
    % for Riemannian submanifolds of the space of matrices), than this
    % function is appropriate:
    % sg = sign(d(find(d(:), 1, 'first')));
	% But it is more difficult to build such a function in general. For
    % now, we ignore this difficulty and let the sign always be +1. This
    % hardly impacts the actual performances, but may be an obstacle for
    % theoretical analysis.
    sg = 1;
    norm_d = problem.M.norm(x, d);
    c = epsilon*sg/norm_d;
    
    % Compute the gradient at the current point.
    [grad0 storedb] = getGradient(problem, x, storedb);
    
    % Compute a point a little further along d and the gradient there.
    x1 = problem.M.retr(x, d, c);
    [grad1 storedb] = getGradient(problem, x1, storedb);
    
    % Transport grad1 back from x1 to x.
    grad1 = problem.M.transp(x1, x, grad1);
    
    % Return the finite difference of them.
    hessfd = problem.M.lincomb(x, 1/c, grad1, -1/c, grad0);
    
end
