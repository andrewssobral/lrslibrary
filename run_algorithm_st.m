%%% ST algorithms
% struct = run_algorithm_st(string, 2dmatrix)
%
function results = run_algorithm_st(algorithm_id, M, opts)
  lrs_load_conf;
  
  alg_path = fullfile(lrs_conf.st_path,algorithm_id);
  addpath(genpath(alg_path));
  
  L = zeros(size(M)); % low-rank matrix
  S = zeros(size(M)); % sparse matrix
  
  results.cputime = 0;
  if(isempty(opts))
    opts.rows = size(M,1);
    opts.cols = size(M,2);
  end
  
  timerVal = tic;
  % warning('off','all');
  
  try
    %
    % ST | GRASTA | Grassmannian Robust Adaptive Subspace Tracking Algorithm (He et al. 2012)
    % process_video('ST', 'GRASTA', 'dataset/demo.avi', 'output/demo_ST-GRASTA.avi');
    %
    if(strcmp(algorithm_id,'GRASTA'))
      run_GRASTA;
    end
    %
    % ST | GOSUS | Grassmannian Online Subspace Updates with Structured-sparsity (Xu et al. 2013)
    % process_video('ST', 'GOSUS', 'dataset/demo.avi', 'output/demo_ST-GOSUS.avi');
    %
    if(strcmp(algorithm_id,'GOSUS'))
      [L,S] = gosus(M,opts);
    end
    %
    % ST | pROST | Robust PCA and subspace tracking from incomplete observations using L0-surrogates (Hage and Kleinsteuber, 2013)
    % process_video('ST', 'pROST', 'dataset/demo.avi', 'output/demo_ST-pROST.avi');
    %
    if(strcmp(algorithm_id,'pROST')) 
      L = robustpca_batch(M,2,'atansquare');
      S = M - L;
    end
  catch ex
    warning(ex.message);
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
