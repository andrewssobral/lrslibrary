%%% NMF algorithms
% struct = run_algorithm_nmf(string, 2dmatrix)
%
function results = run_algorithm_nmf(algorithm_id, M, opts)
  lrs_load_conf;
  
  alg_path = fullfile(lrs_conf.nmf_path,algorithm_id);
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
    % NMF-MU: NMF solved by Multiplicative Updates
    % NMF-PG: NMF solved by Projected Gradient
    % NMF-ALS: NMF solved by Alternating Least Squares
    % NMF-ALS-OBS: NMF solved by Alternating Least Squares with Optimal Brain Surgeon
    % PNMF: Probabilistic Non-negative Matrix Factorization
    %
    % process_video('NMF', 'NMF-MU', 'dataset/demo.avi', 'output/demo_NMF-MU.avi');
    % process_video('NMF', 'NMF-PG', 'dataset/demo.avi', 'output/demo_NMF-PG.avi');
    % process_video('NMF', 'NMF-ALS', 'dataset/demo.avi', 'output/demo_NMF-ALS.avi');
    % process_video('NMF', 'NMF-ALS-OBS', 'dataset/demo.avi', 'output/demo_NMF-ALS-OBS.avi');
    % process_video('NMF', 'PNMF', 'dataset/demo.avi', 'output/demo_PNMF.avi');
    %
    if(strcmp(algorithm_id,'NMF-MU') || strcmp(algorithm_id,'NMF-PG') || ...
      strcmp(algorithm_id,'NMF-ALS') || strcmp(algorithm_id,'NMF-ALS-OBS') || ...
      strcmp(algorithm_id,'PNMF')) 
      alg_path = fullfile(lrs_conf.nmf_path,'NMF-DTU-Toolbox');
      addpath(genpath(alg_path));
      M = sparse(M);
      % mm: Multiplicative update method using euclidean distance measure.
      if(strcmp(algorithm_id,'NMF-MU'))      [W, H] = nmf(M,1,'mm'); end
      % cjlin: Alternative non-negative least squares using projected gradients.
      if(strcmp(algorithm_id,'NMF-PG'))      [W, H] = nmf(M,2,'cjlin'); end
      % als: Alternating least squares.
      if(strcmp(algorithm_id,'NMF-ALS'))     [W, H] = nmf(M,1,'als'); end
      % alsobs: Alternating least squares with optimal brain surgeon.
      if(strcmp(algorithm_id,'NMF-ALS-OBS')) [W, H] = nmf(M,1,'alsobs'); end
      % prob: Probabilistic non-negative matrix factorization.
      if(strcmp(algorithm_id,'PNMF'))    [W, H] = nmf(M,1,'prob'); end
      L = W * H;
      S = M - L;
    end
    %
    % ManhNMF: Manhattan NMF (Guan et al. 2013)
    % process_video('NMF', 'ManhNMF', 'dataset/demo.avi', 'output/demo_ManhNMF.avi');
    if(strcmp(algorithm_id,'ManhNMF'))
      rank = 1;
      [W,H] = ManhNMF(M,rank);
      L = W' * H;
      S = M - L;
    end
    %
    % NeNMF: NMF via Nesterov's Optimal Gradient Method (Guan et al. 2012)
    % process_video('NMF', 'NeNMF', 'dataset/demo.avi', 'output/demo_NeNMF.avi');
    if(strcmp(algorithm_id,'NeNMF'))
      rank = 1;
      [W,H] = NeNMF(M,rank);
      L = W * H;
      S = M - L;
    end
    %
    % LNMF: Spatially Localized NMF (Li et al. 2001)
    % process_video('NMF', 'LNMF', 'dataset/demo.avi', 'output/demo_LNMF.avi');
    if(strcmp(algorithm_id,'LNMF'))
      rank = 1;
      option.verbose = 1;
      [W,H] = LNMF(M,rank,option);
      L = W * H;
      S = M - L;
    end
    %
    % ENMF: Exact NMF (Gillis and Glineur, 2012)
    % process_video('NMF', 'ENMF', 'dataset/demo.avi', 'output/demo_ENMF.avi');
    if(strcmp(algorithm_id,'ENMF'))
      rank = 1;
      [H,W] = ExactNMF(M,rank,100);
      L = (W' * H')';
      S = M - L;
    end
    %
    % nmfLS2: Non-negative Matrix Factorization with sparse matrix (Ji and Eisenstein, 2013)
    % process_video('NMF', 'nmfLS2', 'dataset/demo.avi', 'output/demo_nmfLS2.avi');
    if(strcmp(algorithm_id,'nmfLS2'))
      rank = 1;
      [W,H] = nmfLS2(M, rank);
      L = W * H;
      S = M - L;
    end
    %
    % Semi-NMF: Semi Non-negative Matrix Factorization (Trigeorgis et al. 2014)
    % process_video('NMF', 'Semi-NMF', 'dataset/demo.avi', 'output/demo_Semi-NMF.avi');
    if(strcmp(algorithm_id,'Semi-NMF'))
      rank = 10;
      [W,H] = seminmf(M, rank);
      L = W * H;
      S = M - L;
    end
    %
    % Deep-Semi-NMF: Deep Semi Non-negative Matrix Factorization (Trigeorgis et al. 2014)
    % process_video('NMF', 'Deep-Semi-NMF', 'dataset/demo.avi', 'output/demo_Deep-Semi-NMF.avi');
    if(strcmp(algorithm_id,'Deep-Semi-NMF'))
      layers = 1;
      [W,H] = deep_seminmf(M, layers);
      W = cell2mat(W); H = cell2mat(H);
      L = W * H;
      S = M - L;
    end
    %
    % iNMF: Incremental Subspace Learning via NMF (Bucak and Gunsel, 2009)
    % process_video('NMF', 'iNMF', 'dataset/demo.avi', 'output/demo_iNMF.avi');
    if(strcmp(algorithm_id,'iNMF'))
      % Execute NMF for the first n samples
      n = 10; % first 10 samples
      rdim = 1; % rank-1
      maxiter = 150;
      [W,H] = nmf(M(:,1:n), rdim, 0, maxiter);
      L = W*H;
      % Now we can execute iNMF on each new samples
      maxiter = 50;
      A = M(:,1:n)*H';
      B = H*H';
      h = H(:,end); % Warm start for h
      for i = n+1:size(M,2)
        disp(i);
        M_new = M(:,i);
        [W_new,h,A,B] = inmf(M_new,W,h,A,B,rdim,0.9,0.1,maxiter);
        % H_store(:,i-n) = h; % Just for demonstration
        L(:,end+1) = W_new*h;
      end
      S = M - L;
    end
    %
    % DRMF: Direct Robust Matrix Factorization (Xiong et al. 2011)
    % process_video('NMF', 'DRMF', 'dataset/demo.avi', 'output/demo_DRMF.avi');
    if(strcmp(algorithm_id,'DRMF'))
      %%% initialization
      lambda = 1/sqrt(max(size(M)));
      L_rpca = inexact_alm_rpca(M, lambda, 1e-5, 10);
      sv = svdex(L_rpca);
      rk = EffRank(sv, 0.999);
      %%% run
      options.init = L_rpca;
      [L,S] = DRMF(M, rk, 0.1, options);
      S = full(S);
      %%% end
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
