%function [L,S] = run_SVT(M)
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
I = randperm(numrc)';
k = floor(s*numrc);
I = I(1:k);
MI = M(I);
M2 = M(:); M2(I) = 0;
M2 = reshape(M2,n1,n2);

%%% Set parameters and solve
maxiter = 500;
tol = 1e-4;
tau = 2*5*sqrt(n1*n2); 
delta = 1.2/s;
%{
 if n1 and n2 are very different, then
   tau should probably be bigger than 5*sqrt(n1*n2)

 increase tau to increase accuracy; decrease it for speed

 if the algorithm doesn't work well, try changing tau and delta
   i.e. if it diverges, try a smaller delta (e.g. delta < 2 is a 
   safe choice, but the algorithm may be slower than necessary).
%}

%%% Approximate minimum nuclear norm solution by SVT algorithm
% Note: SVT, as called below, is setup for noiseless data 
%   (i.e. equality constraints).
fprintf('\nSolving by SVT...\n');
%warning('off','SVT:NotUsingMex') 
[U,S,V,numiter] = SVT([n1 n2],I,MI,tau,delta,maxiter,tol);
L = U*S*V';
S = M - L;

% show_2dvideo(M,m,n);
% show_2dvideo(M2,m,n);
% show_2dvideo(L,m,n);
% show_2dvideo(S,m,n);