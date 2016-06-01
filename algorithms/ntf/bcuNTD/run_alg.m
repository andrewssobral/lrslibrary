% NTF | bcuNTD | Non-negative Tucker Decomposition by block-coordinate update (Xu and Yin, 2012)
% process_video('NTF', 'bcuNTD', 'dataset/demo.avi', 'output/demo_bcuNTD.avi');

% Compute a simple NTF model of 10 components
A = double(T);
R = [size(T,1) size(T,2) 2]; % tensor rank
opts.maxit = 1000; % max number of iterations
opts.tol = 1e-4; % stopping tolerance
[M,C] = ntd(T,R,opts);
L = (double(full(ttensor(C,M))));
S = (A - L);
