%% TD algorithms
% struct = run_algorithm_td(string, tensor)
%
function results = run_algorithm_td(algorithm_id, T)
  add_tensor_libs;
  %
  L = []; % low-rank tensor
  S = []; % sparse tensor
  results.cputime = 0;
  %
  % TD | HoSVD | High-order singular value decomposition (Tucker decomposition)
  % process_video('TD', 'HoSVD', 'dataset/demo.avi', 'output/demo_HOSVD.avi');
  %
  if(strcmp(algorithm_id,'HoSVD'))
    % Perform mode-3 rank-1 partial svd
    timerVal = tic;
    [core, U] = tensor_hosvd(T, 0, [0 0 1]);
    L = tensor_ihosvd(core,U);
    cputime = toc(timerVal);
    L = double(L);
    S = double(T) - L;
    % show_3dtensors(T,L,S);
  end
  %
  % TD | t-SVD | Tensor SVD in Fourrier Domain (Zhang et al. 2013)
  % process_video('TD', 't-SVD', 'dataset/demo.avi', 'output/demo_t-SVD.avi');
  %
  if(strcmp(algorithm_id,'t-SVD'))
    A = double(T);
    timerVal = tic;
    [U,S,V] = tensor_t_svd(A);
    [C] = tensor_product(U,S);
    [L] = tensor_product(C,tensor_transpose(V));
    cputime = toc(timerVal);
    S = A - L;
    % show_3dtensors(A,L);
    % show_3dtensors(T,L,S);
  end
  %
  % TD | Tucker-ALS | Tucker Decomposition solved by Alternating Least Squares
  % process_video('TD', 'Tucker-ALS', 'dataset/demo.avi', 'output/demo_Tucker-ALS.avi');
  %
  if(strcmp(algorithm_id,'Tucker-ALS'))
    r = [size(T,1) size(T,2) 1];
    A = double(T);
    timerVal = tic;
    L = double(tucker_als(T,r));
    cputime = toc(timerVal);
    S = A - L;
    % show_3dtensors(A,L);
    % show_3dtensors(T,L,S);
  end
  %
  % TD | CP-ALS | PARAFAC/CP decomposition solved by Alternating Least Squares
  % process_video('TD', 'CP-ALS', 'dataset/demo.avi', 'output/demo_CP-ALS.avi');
  %
  if(strcmp(algorithm_id,'CP-ALS'))
    r = 10;
    A = double(T);
    timerVal = tic;
    L = double(cp_als(T,r,'dimorder',[3 2 1]));
    cputime = toc(timerVal);
    S = A - L;
    % show_3dtensors(A,L);
    % show_3dtensors(T,L,S);
  end
  %
  % TD | CP-APR | PARAFAC/CP decomposition solved by Alternating Poisson Regression
  % process_video('TD', 'CP-APR', 'dataset/demo.avi', 'output/demo_CP-APR.avi');
  %
  if(strcmp(algorithm_id,'CP-APR'))
    r = 10;
    A = double(T);
    timerVal = tic;
    L = double(cp_apr(T,r));
    cputime = toc(timerVal);
    S = A - L;
    % show_3dtensors(A,L);
    % show_3dtensors(T,L,S);
  end
  %
  % TD | CP2 | PARAFAC2 decomposition
  % process_video('TD', 'CP2', 'dataset/demo.avi', 'output/demo_CP2.avi');
  %
  if(strcmp(algorithm_id,'CP2'))
    addpath('algorithms/td/parafac2');
    r = 10;
    A = double(T);
    timerVal = tic;
    [B,H,C,P] = parafac2(A,r,[],[0 0 0 0 1]);
    %%% PARAFAC2 reconstruction
    for i = 1:size(C,1)
      L(:,:,i) = B*diag(C(i,:))*(P{i}*H)';
    end
    cputime = toc(timerVal);
    S = A - L;
    % show_3dtensors(A,L);
    % show_3dtensors(T,L,S);
    rmpath('algorithms/td/parafac2');
  end
  %
  % TD | RSTD | Rank Sparsity Tensor Decomposition (Yin Li 2010)
  % process_video('TD', 'RSTD', 'dataset/demo.avi', 'output/demo_RSTD.avi');
  %
  if(strcmp(algorithm_id,'RSTD'))
    addpath('algorithms/td/rstd');
    addpath('algorithms/td/rstd/utils');
    addpath('libs/PROPACK');
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
    timerVal = tic;
    [L, S] = RSTD(A, alpha, beta, gamma, lambda, eta, maxIter, rank);
    cputime = toc(timerVal);
    % show_3dtensors(A,L);
    % show_3dtensors(T,L,S);
    rmpath('libs/PROPACK');
    rmpath('algorithms/td/rstd/utils');
    rmpath('algorithms/td/rstd');
  end
  %
  % TD | HoRPCA-IALM | HoRPCA solved by IALM (Goldfarb and Qin, 2013)
  % TD | HoRPCA-S | HoRPCA with Singleton model solved by ADAL (Goldfarb and Qin, 2013)
  % TD | HoRPCA-S-NCX | HoRPCA with Singleton model solved by ADAL (non-convex) (Goldfarb and Qin, 2013)
  % TD | Tucker-ADAL | Tucker Decomposition solved by ADAL (Goldfarb and Qin, 2013)
  %
  % process_video('TD', 'HoRPCA-IALM', 'dataset/demo.avi', 'output/demo_HoRPCA-IALM.avi');
  % process_video('TD', 'HoRPCA-S', 'dataset/demo.avi', 'output/demo_HoRPCA-S.avi');
  % process_video('TD', 'Tucker-ADAL', 'dataset/demo.avi', 'output/demo_Tucker-ADAL.avi');
  %
  if(strcmp(algorithm_id,'HoRPCA-IALM') || strcmp(algorithm_id,'HoRPCA-S') ...
  || strcmp(algorithm_id,'HoRPCA-S-NCX')|| strcmp(algorithm_id,'Tucker-ADAL'))
    addpath(genpath('algorithms/td/rlrt'));
    addpath('libs/PROPACK');
    data.T = T;
    data.X = T;
    N = ndims(data.T);
    r = 1/sqrt(max(size(data.T))); 
    params.E0 = tenzeros(size(data.T));
    params.X0 = tenzeros(size(data.T));
    params.V0 = cell(1, N);
    for i = 1:N 
      params.V0{i} = tenzeros(size(data.T));
    end
    params.mu0 = 1/(N+1);
    params.mode = N;
    params.IsTC = false; % is tensor completion
    params.rRatio = 1/4;
    params.opt_tol = 1e-3;
    params.eta = 1/(N+1);
    params.max_iter = 1000;
    params.mu1fac = 10;
    params.mu1 = params.mu1fac*std(T(:));
    params.mu2 = params.mu1;
    params.mu_min = 1e-4;
    params.mu_max = 1e2;
    params.lambdaS = 1;
    params.lambda = params.lambdaS*r*params.rRatio; 
    params.verbose = 1;
    params.use_cont = true;
    params.k = [size(T,1) size(T,2) 1];
    %%%%%%%%%% for PROPACK %%%%%%%%%%%%
    % declare global var 'sv'
    global sv;
    global tmode;
    global use_propack;
    global curr_mu;
    sv =  ceil(min(size(data.T)) * 0.1) * ones( 1, N );
    use_propack = true;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    timerVal = tic;
    if(strcmp(algorithm_id,'HoRPCA-IALM')) results = rpca_for_tensor(data, params); end
    if(strcmp(algorithm_id,'HoRPCA-S')) results = tensor_rpca_adal2(data, params); end % tensor_rpca_adal
    if(strcmp(algorithm_id,'HoRPCA-S-NCX')) results = tensor_rpca_adal_ncx(data, params); end
    if(strcmp(algorithm_id,'Tucker-ADAL')) results = tensor_tucker_adal_ncx(data, params); end
    cputime = toc(timerVal);
    L = double(results.X);
    S = double(results.E);
    % show_3dtensors(A,L);
    % show_3dtensors(T,L,S);
    clear sv tmode use_propack curr_mu params results data N;
    rmpath('libs/PROPACK');
    rmpath(genpath('algorithms/td/rlrt'));
  end
  %
  %
  %
  results.L = L; % low-rank matrix
  results.S = S; % sparse matrix
  results.cputime = cputime;
  clear L S;
end
