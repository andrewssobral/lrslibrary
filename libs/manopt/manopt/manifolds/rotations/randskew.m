function S = randskew(n, N)
% Generates random skew symmetric matrices with normal entries.
% 
% function S = randskew(n, N)
%
% S is an n-by-n-by-N matrix where each slice S(:, :, i) for i = 1..N is a
% random skew-symmetric matrix with upper triangular entries distributed
% independently following a normal distribution (Gaussian, zero mean, unit
% variance).
%
% See also: randrot

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Sept. 25, 2012.
% Contributors: 
% Change log: 


    if nargin < 2
        N = 1;
    end

    % Subindices of the (strictly) upper triangular entries of an n-by-n
    % matrix
    [I J] = find(triu(ones(n), 1));
    
    K = repmat(1:N, n*(n-1)/2, 1);
    
    % Indices of the strictly upper triangular entries of all N slices of
    % an n-by-n-by-N matrix
    L = sub2ind([n n N], repmat(I, N, 1), repmat(J, N, 1), K(:));
    
    % Allocate memory for N random skew matrices of size n-by-n and
    % populate each upper triangular entry with a random number following a
    % normal distribution and copy them with opposite sign on the
    % corresponding lower triangular side.
    S = zeros(n, n, N);
    S(L) = randn(size(L));
    S = S-multitransp(S);

end
