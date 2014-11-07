%%% NMF algorithms
% struct = run_algorithm_nmf(string, 2dmatrix)
%
function results = run_algorithm_nmf(algorithm_id, M, opts)
  L = []; % low-rank matrix
  S = []; % sparse matrix
  results.cputime = 0;
  if(isempty(opts))
    opts.rows = size(M,1);
    opts.cols = size(M,2);
  end
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
    addpath('algorithms/nmf/NMF-DTU-Toolbox');
    M = sparse(M);
    timerVal = tic;
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
    cputime = toc(timerVal);
    L = W * H;
    S = M - L;
    % imagesc(L); imagesc(S);
    clear W H;
    rmpath('algorithms/nmf/NMF-DTU-Toolbox');
  end
  %
  % ManhNMF: Manhattan NMF (Guan et al. 2013)
  % process_video('NMF', 'ManhNMF', 'dataset/demo.avi', 'output/demo_ManhNMF.avi');
  if(strcmp(algorithm_id,'ManhNMF')) 
    addpath('algorithms/nmf/ManhNMF');
    rank = 1;
    timerVal = tic;
    [W, H] = ManhNMF(M,rank);
    cputime = toc(timerVal);
    L = W' * H;
    S = M - L;
    % imagesc(L); imagesc(S);
    clear rank W H;
    rmpath('algorithms/nmf/ManhNMF');
  end
  %
  % NeNMF: NMF via Nesterov's Optimal Gradient Method (Guan et al. 2012)
  % process_video('NMF', 'NeNMF', 'dataset/demo.avi', 'output/demo_NeNMF.avi');
  if(strcmp(algorithm_id,'NeNMF')) 
    addpath('algorithms/nmf/NeNMF');
    rank = 1;
    timerVal = tic;
    [W, H] = NeNMF(M,rank);
    cputime = toc(timerVal);
    L = W * H;
    S = M - L;
    % imagesc(L); imagesc(S);
    clear rank W H;
    rmpath('algorithms/nmf/NeNMF');
  end
  %
  % LNMF: Spatially Localized NMF (Li et al. 2001)
  % process_video('NMF', 'LNMF', 'dataset/demo.avi', 'output/demo_LNMF.avi');
  if(strcmp(algorithm_id,'LNMF')) 
    addpath('algorithms/nmf/LNMF');
    rank = 1;
    option.verbose = 1;
    timerVal = tic;
    [W, H] = LNMF(M,rank,option);
    cputime = toc(timerVal);
    L = W * H;
    S = M - L;
    % imagesc(L); imagesc(S);
    clear rank W H;
    rmpath('algorithms/nmf/LNMF');
  end
  %
  % ENMF: Exact NMF (Gillis and Glineur, 2012)
  % process_video('NMF', 'ENMF', 'dataset/demo.avi', 'output/demo_ENMF.avi');
  if(strcmp(algorithm_id,'ENMF')) 
    addpath('algorithms/nmf/ENMF');
    rank = 1;
    timerVal = tic;
    [H, W] = ExactNMF(M,rank,100);
    cputime = toc(timerVal);
    L = (W' * H')';
    S = M - L;
    % imagesc(L); imagesc(S);
    clear rank W H;
    rmpath('algorithms/nmf/ENMF');
  end
  %
  % nmfLS2: Non-negative Matrix Factorization with sparse matrix (Ji and Eisenstein, 2013)
  % process_video('NMF', 'nmfLS2', 'dataset/demo.avi', 'output/demo_nmfLS2.avi');
  if(strcmp(algorithm_id,'nmfLS2')) 
    addpath('algorithms/nmf/nmfLS2');
    rank = 1;
    timerVal = tic;
    [W, H] = nmfLS2(M, rank);
    cputime = toc(timerVal);
    L = W * H;
    S = M - L;
    % imagesc(L); imagesc(S);
    clear rank W H;
    rmpath('algorithms/nmf/nmfLS2');
  end
  %
  % Semi-NMF: Semi Non-negative Matrix Factorization (Trigeorgis et al. 2014)
  % process_video('NMF', 'Semi-NMF', 'dataset/demo.avi', 'output/demo_Semi-NMF.avi');
  if(strcmp(algorithm_id,'Semi-NMF')) 
    addpath('algorithms/nmf/Semi-NMF');
    rank = 10;
    timerVal = tic;
    [W, H] = seminmf(M, rank);
    cputime = toc(timerVal);
    L = W * H;
    S = M - L;
    % clf; imagesc(L); imagesc(S);
    clear rank W H;
    rmpath('algorithms/nmf/Semi-NMF');
  end
  %
  % Deep-Semi-NMF: Deep Semi Non-negative Matrix Factorization (Trigeorgis et al. 2014)
  % process_video('NMF', 'Deep-Semi-NMF', 'dataset/demo.avi', 'output/demo_Deep-Semi-NMF.avi');
  if(strcmp(algorithm_id,'Deep-Semi-NMF')) 
    addpath('algorithms/nmf/Deep-Semi-NMF');
    layers = 1;
    timerVal = tic;
    [W, H] = deep_seminmf(M, layers);
    cputime = toc(timerVal);
    W = cell2mat(W); H = cell2mat(H);
    L = W * H;
    S = M - L;
    % clf; imagesc(L); imagesc(S);
    clear rank W H;
    rmpath('algorithms/nmf/Deep-Semi-NMF');
  end
  %
  %
  %
  results.L = L; % low-rank matrix
  results.S = S; % sparse matrix
  results.cputime = cputime;
  clear L S;
end
