%%% TTD algorithms
% struct = run_algorithm_ttd(string, 2dmatrix)
%
function results = run_algorithm_ttd(algorithm_id, M, opts)
  lrs_load_conf;
  
  alg_path = fullfile(lrs_conf.ttd_path,algorithm_id);
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
    % TTD | 3WD | 3-Way-Decomposition (Oreifej et al. 2012)
    % process_video('TTD', '3WD', 'dataset/demo.avi', 'output/demo_3WD.avi');
    %
    if(strcmp(algorithm_id,'3WD'))
      run_3WD;
    end
    %
    % TTD | MAMR | Motion-Assisted Matrix Restoration (Ye et al. 2015)
    % process_video('TTD', 'MAMR', 'dataset/demo.avi', 'output/demo_MAMR.avi');
    %
    if(strcmp(algorithm_id,'MAMR') || strcmp(algorithm_id,'RMAMR'))
      alg_path = fullfile(lrs_conf.ttd_path,'MAMR_RMAMR');
      addpath(genpath(alg_path));
      run_MAMR_RMAMR;
    end
    %
    % TTD | ADMM | Alternating Direction Method of Multipliers (Parikh and Boyd, 2014)
    % process_video('TTD', 'ADMM', 'dataset/demo.avi', 'output/demo_ADMM.avi');
    %
    if(strcmp(algorithm_id,'ADMM'))
      results = ADMM(M); % show_2dvideo(M,m,n);
      %E = results.Z; % show_2dvideo(E,m,n);
      S = results.S; % show_2dvideo(S,m,n);
      L = results.L; % show_2dvideo(L,m,n);
      %M_hat = L + S + E; % show_2dvideo(M_hat,m,n);
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
