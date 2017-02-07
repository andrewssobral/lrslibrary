%{
load('dataset/trafficdb/traffic_patches.mat');
[M,m,n,p] = convert_video3d_to_2d(im2double(imgdb{100}));
out = run_algorithm('MC', 'SVT', M, [])
show_results(M.*out.Omega,out.L,out.S,out.O,p,m,n);
%}

MIdx = M(Idx);
maxiter = 500;
tol = 1e-4;
tau = 2*5*sqrt(numel(M));
delta = 1.2; %/s;

%%% Approximate minimum nuclear norm solution by SVT algorithm
% Note: SVT, as called below, is setup for noiseless data
%   (i.e. equality constraints).
fprintf('\nSolving by SVT...\n');
%warning('off','SVT:NotUsingMex')
[U,S,V,numiter] = SVT(size(M),Idx,MIdx,tau,delta,maxiter,tol);
L = U*S*V'; % low-rank
S = M - L; % sparse
