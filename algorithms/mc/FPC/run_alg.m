% MC | FPC | Fixed point and Bregman iterative methods for matrix rank minimization (Ma et al. 2008)
% process_video('MC', 'FPC', 'dataset/demo.avi', 'output/demo_MC-FPC.avi');
%{
load('dataset/trafficdb/traffic_patches.mat');
[M,m,n] = convert_video3d_to_2d(im2double(imgdb{100}));
%}

%sigma = .05*std(M(:));
%M = M + sigma*randn(size(M));

%M = (M - mean(M(:)))./(var(M(:)));
%M(M >= 0.999) = 0.999; % best
%M(M <= 0.001) = 0.001; % best
%M = M./norm(M);
%M = M./norm(M,1);
%M = M/max(abs(M(:)));
%min(M(:))
%max(M(:))
%M = nma_rescale(M,0.207,0.968);
M = nma_rescale(M,0.01,0.99);

s = 0.5; % 50% of observed values
[n1,n2] = size(M);
numrc = n1*n2;
Idx = randperm(numrc)';
k = floor(s*numrc);
Idx = Idx(1:k);
MI = M(Idx);
M2 = M(:); M2(Idx) = 0;
M2 = reshape(M2,n1,n2);

% Approximate minimum nuclear norm solution by FPC algorithm
% This version of FPC uses PROPACK for the multiplies
%   It is not optimized, and the parameters have not been tested
% The version in the FPC paper by Shiqian Ma, Donald Goldfarb and Lifeng Chen
%   uses an approximate SVD that will have different properties;
%   That code may be found at http://www.columbia.edu/~sm2756/FPCA.htm
maxiter = 500;
mu_final = .01; tol = 1e-3;

fprintf('\nSolving by FPC...\n');
[U,S,V,numiter] = FPC([n1 n2],Idx,MI,mu_final,maxiter,tol);

L = U*S*V';
S = M - L;

% show_2dvideo(M,m,n);
% show_2dvideo(M2,m,n);
% show_2dvideo(L,m,n);
% show_2dvideo(S,m,n);