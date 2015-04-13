function R = randrot(n, N)
% Generates uniformly random rotation matrices.
%
% function R = randrot(n, N)
%
% R is a n-by-n-by-N matrix such that each slice R(:, :, i) is an
% orthogonal matrix of size n of determinant +1 (i.e., a matrix in SO(n)).
% By default, N = 1.
% Complexity: N times O(n^3).
% Theory in Diaconis and Shahshahani 1987 for the uniformity on O(n);
% With details in Mezzadri 2007,
% "How to generate random matrices from the classical compact groups."
% To ensure matrices in SO(n), we permute the two first columns when
% the determinant is -1.
%
% See also: randskew

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Sept. 25, 2012.
% Contributors: 
% Change log: 

    if nargin < 2
        N = 1;
    end
    
    if n == 1
        R = ones(1, 1, N);
        return;
    end
    
    R = zeros(n, n, N);
    
    for i = 1 : N
        
        % Generated as such, Q is uniformly distributed over O(n), the set
        % of orthogonal matrices.
        A = randn(n);
        [Q, RR] = qr(A);
        Q = Q * diag(sign(diag(RR))); %% Mezzadri 2007
        
        % If Q is in O(n) but not in SO(n), we permute the two first
        % columns of Q such that det(new Q) = -det(Q), hence the new Q will
        % be in SO(n), uniformly distributed.
        if det(Q) < 0
            Q(:, [1 2]) = Q(:, [2 1]);
        end
        
        R(:, :, i) = Q;
        
    end

end
