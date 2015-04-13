%%% MC algorithms
% struct = run_algorithm_mc(string, 2dmatrix)
%
function results = run_algorithm_mc(algorithm_id, M, opts)
  lrs_load_conf;
  
  alg_path = fullfile(lrs_conf.mc_path,algorithm_id);
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
  
  %try
    %
    % MC | LRGeomCG | Low-rank matrix completion by Riemannian optimization (Bart Vandereycken, 2013)
    % process_video('MC', 'LRGeomCG', 'dataset/demo.avi', 'output/demo_MC-LRGeomCG.avi');
    %
    if(strcmp(algorithm_id,'LRGeomCG'))
      [L,S] = low_rank_matrix_completion(M);
    end
    %
    % MC | GROUSE | Grassmannian Rank-One Update Subspace Estimation (Balzano et al. 2010)
    % process_video('MC', 'GROUSE', 'dataset/demo.avi', 'output/demo_MC-GROUSE.avi');
    %
    if(strcmp(algorithm_id,'GROUSE'))
      run_GROUSE;
    end
    %
    % MC | OptSpace | A Matrix Completion Algorithm (Keshavan et al. 2009)
    % process_video('MC', 'OptSpace', 'dataset/demo.avi', 'output/demo_MC-OptSpace.avi');
    %
    if(strcmp(algorithm_id,'OptSpace'))
      run_OptSpace;
    end
    %
    % MC | FPC | Fixed point and Bregman iterative methods for matrix rank minimization (Ma et al. 2008)
    % process_video('MC', 'FPC', 'dataset/demo.avi', 'output/demo_MC-FPC.avi');
    %
    if(strcmp(algorithm_id,'FPC'))
      run_FPC;
    end
    %
    % MC | SVT | A singular value thresholding algorithm for matrix completion (Cai et al. 2008)
    % process_video('MC', 'SVT', 'dataset/demo.avi', 'output/demo_MC-SVT.avi');
    %
    if(strcmp(algorithm_id,'SVT'))
      run_SVT;
    end
  %catch ex
  %  warning(ex.message);
  %end
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
