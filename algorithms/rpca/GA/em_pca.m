% EM_PCA    Estimate principal components using an EM algorithm.
%
%     EM_PCA(X) returns a basis vector for the first principal subspace of X
%     This is a N-by-D matrix containing N observations in R^D.
%
%     EM_PCA(X, K) returns a D-by-K matrix of orthogonal principal vectors.
%
%     Reference:
%     S.T. Roweis. EM algorithms for PCA and SPCA. NIPS, 1998.

%    Grassmann Averages
%    Copyright (C) 2014  Søren Hauberg
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function vectors = em_pca(X, K)
  %% Check input
  if (nargin == 0)
    error('em_pca: not enough input arguments');
  elseif (nargin == 1)
    K = 1;
  end % if
  
  epsilon = 1e-4; % Note that smaller values are highly impractical
  
  %% Create output structure
  [N, D] = size(X);
  vectors = NaN(D, K, class(X));
  
  for k = 1:K
    %% Compute k'th principal component
    mu = rand(D, 1) - 0.5;
    mu = mu / norm(mu);
    for iter = 1:N
      prev_mu = mu;

      %% EM step (note similarity to Grassmann Average)
      dots = X * mu;
      mu = dots.' * X;
      mu = mu(:) / norm(mu);

      %% Check for convergence
      if (sum(abs(mu - prev_mu)) < epsilon)
        break;
      end % if
    end % for

    %% Store the estimated vector (and possibly subtract it from data, and perform reorthonomralisation)
    if (k == 1)
      vectors(:, k) = mu;
    elseif (k < K)
      mu = reorth(vectors(:, 1:k-1), mu, 1);
      mu = mu / norm(mu);
      vectors(:, k) = mu;

      X = X - (X * mu) * mu.';
    else % k == K
      mu = reorth(vectors(:, 1:k-1), mu, 1);
      mu = mu / norm(mu);
      vectors(:, k) = mu;
    end % if
  end % for
end % function

