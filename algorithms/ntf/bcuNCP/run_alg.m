% NTF | bcuNCP | Non-negative CP Decomposition by block-coordinate update (Xu and Yin, 2012)
% process_video('NTF', 'bcuNCP', 'dataset/demo.avi', 'output/bcuNCP.avi');

% Compute a simple NTF model of 10 components
A = double(T);
R = 10; % tensor rank
opts.maxit = 1000; % max number of iterations
opts.tol = 1e-4; % stopping tolerance
M = ncp(T,R,opts);
L = double(full(M));
S = (A - L);
