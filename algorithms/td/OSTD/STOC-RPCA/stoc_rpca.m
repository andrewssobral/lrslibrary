% Stochastic optimization for the robust PCA
% Input:
%       D: [m x n] data matrix, m - ambient dimension, n - samples number
%       lambda1, lambda2: trade-off parameters
%       nrank: the groundtruth rank of the data
% Output:
%       L: [m x r] the basis of the subspace
%       R: [r x n] the coefficient of the samples on the basis L
%       E: sparse error
%
% copyright (c) Jiashi Feng (jshfeng@gmail.com)
%
% This is a modified version of Jiashi Feng code
%
%{
lambda1 = 1/sqrt(max(size(M)));
lambda2 = lambda1;
nrank = size(M,2);
[L,S] = stoc_rpca(M, lambda1, lambda2, nrank);
%}
% [L,R,E] = stoc_rpca(D, lambda1, lambda2, nrank)
function [Lr,S] = stoc_rpca(D, lambda1, lambda2, nrank)
  %% initialization
  [ndim, nsample] = size(D);
  %L = cell(nsample+1,1);
  %L{1} = rand(ndim,nrank);
  L = rand(ndim,nrank);
  A = zeros(nrank,nrank);
  B = zeros(ndim,nrank);
  %R = zeros(nrank,nsample);
  S = zeros(ndim,nsample);
  Lr = zeros(ndim,nsample);
  %% online optimization
  for t = 1:nsample
    disp(t);
    %   tic;
    z = D(:,t);
    %   tic;
    %[r,s] = solve_proj2(z,L{t},lambda1,lambda2);
    [r,s] = solve_proj2(z,L,lambda1,lambda2);
    %   tused = toc;
    %   fprintf('elapsed time for projection %f secs\n',tused);
    %   R(:,t) = r;
    S(:,t) = s;
    A = A + r*r';
    B = B + (z-s)*r';
    %   L{t+1} = update_col(L{t},A,B,lambda1/nsample);
    L = update_col_orpca(L,A,B,lambda1/nsample);
    Lr(:,t) = L * r;
    %   disp(t);
    %   toc;
  end
end
