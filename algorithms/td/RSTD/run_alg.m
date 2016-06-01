% TD | RSTD | Rank Sparsity Tensor Decomposition (Yin Li 2010)
% process_video('TD', 'RSTD', 'dataset/demo.avi', 'output/demo_RSTD.avi');

maxIter = 400;                        % maximun iteration number
alpha = [1, 1, 0.1];                  % relaxation parameter for rank
beta = [1, 1, 0.1];                   % relaxation parameter for sparsity
gamma = [1, 1, 0.1];                  % relaxation parameter for consistency
lambda = [4.8, 4.8, 0.1];             % the weights of trace norm terms
eta = [0.1, 0.1, 0.1];                % the weights of l1 norm terms
A = double(T);
rank = [size(A,1) size(A,2) 1];
%[TL, TS, Ud, rank2, sparsity, errorList, iter] = ...
%  RSTD(A, alpha, beta, gamma, lambda, eta, maxIter, rank);
%core = HOSVD(A, Ud);
%A_hat = iHOSVD(core, Ud);
[L,S] = RSTD(A, alpha, beta, gamma, lambda, eta, maxIter, rank);
