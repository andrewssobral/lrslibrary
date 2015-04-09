%%% LRR algorithms
% struct = run_algorithm_lrr(string, 2dmatrix)
%
function results = run_algorithm_lrr(algorithm_id, M, opts)
  load_lrslibrary_conf;
  
  alg_path = fullfile(lrs_conf.lrr_path,algorithm_id);
  addpath(genpath(alg_path));
  
  L = []; % low-rank matrix
  S = []; % sparse matrix
  
  results.cputime = 0;
  if(isempty(opts))
    opts.rows = size(M,1);
    opts.cols = size(M,2);
  end
  
  timerVal = tic;
  % warning('off','all');
  
  %
  % LRR | EALM | Exact ALM (Lin et al. 2009)
  % LRR | IALM | Inexact ALM (Lin et al. 2009)
  %
  % process_video('LRR', 'EALM', 'dataset/demo.avi', 'output/demo_LRR-EALM.avi');
  % process_video('LRR', 'IALM', 'dataset/demo.avi', 'output/demo_LRR-IALM.avi');
  %
  if(strcmp(algorithm_id,'EALM') || strcmp(algorithm_id,'IALM'))
    alg_path = fullfile(lrs_conf.lrr_path,'ALM');
    addpath(genpath(alg_path));
    A = mean(M,2);
    lambda = 0.01;
    if(strcmp(algorithm_id,'EALM')) [Z,E] = solve_lrr(M,A,lambda,0,0,1); end
    if(strcmp(algorithm_id,'IALM')) [Z,E] = solve_lrr(M,A,lambda,0,1,1); end
    % M_hat = A*Z + E;
    L = A*Z;
    S = E;
  end
  %
  % LRR | ADM | Alternating Direction Method (Lin et al. 2011)
  % LRR | LADMAP | Linearized ADM with Adaptive Penalty (Lin et al. 2011)
  % LRR | FastLADMAP | Fast LADMAP (Lin et al. 2011)
  %
  % process_video('LRR', 'ADM', 'dataset/demo.avi', 'output/demo_LRR-ADM.avi');
  % process_video('LRR', 'LADMAP', 'dataset/demo.avi', 'output/demo_LRR-LADMAP.avi');
  % process_video('LRR', 'FastLADMAP', 'dataset/demo.avi', 'output/demo_LRR-FastLADMAP.avi');
  %
  if(strcmp(algorithm_id,'ADM') || strcmp(algorithm_id,'LADMAP') || strcmp(algorithm_id,'FastLADMAP'))
    alg_path = fullfile(lrs_conf.lrr_path,'ADM');
    addpath(genpath(alg_path));
    lambda = 0.1;
    rho = 1.9;
    DEBUG = 1;
    % X = XZ+E
    if(strcmp(algorithm_id,'ADM'))        [Z,E] = adm_lrr(M,lambda,rho,DEBUG); end
    if(strcmp(algorithm_id,'LADMAP'))     [Z,E] = ladmp_lrr(M,lambda,rho,DEBUG); clearvars -global M; end
    if(strcmp(algorithm_id,'FastLADMAP')) [Z,E] = ladmp_lrr_fast(M,lambda,rho,DEBUG); clearvars -global A Xg eta M; end
    %M_hat = M*Z + E;
    L = M*Z;
    S = E; %S = M_hat - L;
  end
  %
  % LRR | ROSL | Robust Orthonormal Subspace Learning (Shu et al. 2014)
  % process_video('LRR', 'ROSL', 'dataset/demo.avi', 'output/demo_LRR-ROSL.avi');
  %
  if(strcmp(algorithm_id,'ROSL'))
    K = 1; % The initialiation of the subspace dimension
    tol = 1e-5;
    maxIter = 30;
    lambda = 1e-1; %2e-3;
    [~,~,E_hat,A_hat] = inexact_alm_rosl(M,K,lambda,tol,maxIter);
    L = A_hat;
    S = E_hat;
  end
  %
  cputime = toc(timerVal);
  rmpath(genpath(alg_path));
  %
  results.L = L; % low-rank matrix
  results.S = S; % sparse matrix
  results.cputime = cputime;
  %
  % warning('on','all');
end
