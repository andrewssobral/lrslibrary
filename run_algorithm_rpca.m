%%% RPCA algorithms
% struct = run_algorithm_rpca(string, 2dmatrix)
%
function results = run_algorithm_rpca(algorithm_id, M, opts)
  L = []; % low-rank matrix
  S = []; % sparse matrix
  results.cputime = 0;
  if(isempty(opts))
    opts.rows = size(M,1);
    opts.cols = size(M,2);
  end
  %
  % PCP (Candes et al. 2009)
  % process_video('RPCA', 'PCP', 'dataset/demo.avi', 'output/demo_PCP.avi');
  if(strcmp(algorithm_id,'PCP')) 
    addpath('algorithms/rpca/PCP');
    lambda = 1/sqrt(max(size(M,1),size(M,2)));
    tol = 1e-5;
    timerVal = tic;
    [L, S] = PCP(M,lambda,tol); % imagesc(L); imagesc(S);
    cputime = toc(timerVal);
    clear lambda tol;
    rmpath('algorithms/rpca/PCP');
  end
  %
  % Fast PCP (Rodriguez and Wohlberg 2013)
  % process_video('RPCA', 'FPCP', 'dataset/demo.avi', 'output/demo_FPCP.avi');
  if(strcmp(algorithm_id,'FPCP')) 
    addpath('algorithms/rpca/FastPCP');
    lambda = 1/sqrt(max(size(M,1),size(M,2)));
    timerVal = tic;
    [L, S] = fastpcp(M, lambda); % imagesc(L); imagesc(S);
    cputime = toc(timerVal);
    clear lambda;
    rmpath('algorithms/rpca/FastPCP');
  end
  %
  % AS-RPCA: Active Subspace: Towards Scalable Low-Rank Learning (Liu and Yan, 2012)
  % process_video('RPCA', 'AS-RPCA', 'dataset/demo.avi', 'output/demo_AS-RPCA.avi');
  if(strcmp(algorithm_id,'AS-RPCA')) 
    addpath('algorithms/rpca/AS-RPCA');
    timerVal = tic;
    lambda = 1/sqrt(min(size(M,1),size(M,2)));
    [L, S] = as_rpca(M,lambda); % imagesc(L); imagesc(S);
    cputime = toc(timerVal);
    clear lambda;
    rmpath('algorithms/rpca/AS-RPCA');
  end
  %
  % R2PCP: Riemannian Robust Principal Component Pursuit (Hintermüller and Wu, 2014)
  % process_video('RPCA', 'R2PCP', 'dataset/demo.avi', 'output/demo_R2PCP.avi');
  if(strcmp(algorithm_id,'R2PCP')) 
    addpath('algorithms/rpca/R2PCP');
    addpath('libs/PROPACK');
    load_sig = 0; % loading signal
    input_sig = 2;
    noise = 0e-4; % noise level
    % data_formation
    if ~load_sig, A_ref=[]; B_ref=[]; end
    Z = M;
    %[z1,z2,z3] = size(Z);
    z1 = opts.rows;
    z2 = opts.cols;
    z3 = size(Z,2);
    n1 = z1*z2; n2 = z3;
    %%% model setting
    model_formation;
    %%% algorithm parameters
    slv.kmax=100;           % AMS max iteration
    slv.ktol=1e-4;          % tolerance on residual norm
    ls.imax=20;             % line search max iteration
    % ls.c=.001;
    cg.imax=30;             % CG max iteration
    cg.tol=1e-2;            % CG residual tolerance
    trim.sig=2;             % trimming signal (=2 recommended)
    trim.disp=0;
    trim.tol=2e-1;
    trim.mod=1;             % frequency of trimming
    trim.thr1=log10(10);    % threshold for singular values of A
    trim.thr2=.12;          % threshold for absolute values of B
    sg.sig=1;               % safeguard signal:
                            % 0: no safteguard
                            % 1: safeguard for B-subprob only (A-subprob via projected dogleg)
                            % 2: safeguard for both subproblems
    sg.thr=1e-1;            % threshold for safeguard
    riop.sig=2;             % 1: Riemannian gradient
                            % 2: projected dogleg
    riop.eps_c=0e-6;
    %%% initialization
    n3 = 5;                   % initial rank of A
    n4 = round(n1*n2*.1);     % initial cardinality of B
    initialization;
    % alternating minimization scheme
    timerVal = tic;
    slv_lrs_ams;
    cputime = toc(timerVal);
    L = A; % imagesc(L);
    S = B; % imagesc(S);
    rmpath('libs/PROPACK');
    rmpath('algorithms/rpca/R2PCP');
  end
  %
  % ALM (Tang and Nehorai 2011)
  % process_video('RPCA', 'ALM', 'dataset/demo.avi', 'output/demo_ALM.avi');
  if(strcmp(algorithm_id,'ALM')) 
    addpath('algorithms/rpca/ALM');
    addpath('libs/PROPACK');
    timerVal = tic;
    [L, S] = alm(M); % imagesc(L); imagesc(S);
    cputime = toc(timerVal);
    rmpath('libs/PROPACK');
    rmpath('algorithms/rpca/ALM');
  end
  %
  % EALM (Lin et al. 2009)
  % process_video('RPCA', 'EALM', 'dataset/demo.avi', 'output/demo_EALM.avi');
  if(strcmp(algorithm_id,'EALM')) 
    addpath('algorithms/rpca/EALM');
    timerVal = tic;
    [L, S] = exact_alm_rpca(M); % imagesc(L); imagesc(S);
    cputime = toc(timerVal);
    rmpath('algorithms/rpca/EALM');
  end
  %
  % IALM (Lin et al. 2009)
  % process_video('RPCA', 'IALM', 'dataset/demo.avi', 'output/demo_IALM.avi');
  if(strcmp(algorithm_id,'IALM')) 
    addpath('algorithms/rpca/IALM');
    timerVal = tic;
    [L, S] = inexact_alm_rpca(M); % imagesc(L); imagesc(S);
    cputime = toc(timerVal);
    rmpath('algorithms/rpca/IALM');
  end
  %
  % IALM + LMSVDS (Liu et al. 2012)
  % process_video('RPCA', 'IALM_LMSVDS', 'dataset/demo.avi', 'output/demo_IALM_LMSVDS.avi');
  if(strcmp(algorithm_id,'IALM_LMSVDS'))
    warning('off','all');
    addpath('algorithms/rpca/IALM');
    timerVal = tic;
    [L, S] = inexact_alm_rpca_with_lmsvds(M); % imagesc(L); imagesc(S);
    cputime = toc(timerVal);
    rmpath('algorithms/rpca/IALM');
    warning('on','all');
  end
  %
  % IALM + BLWS (Lin and Wei 2010)
  % process_video('RPCA', 'IALM_BLWS', 'dataset/demo.avi', 'output/demo_IALM_BLWS.avi');
  if(strcmp(algorithm_id,'IALM_BLWS'))
    warning('off','all');
    addpath('algorithms/rpca/IALM');
    addpath('libs/PROPACK');
    timerVal = tic;
    [L, S] = inexact_alm_rpca_with_blws(M); % imagesc(L); imagesc(S);
    cputime = toc(timerVal);
    rmpath('libs/PROPACK');
    rmpath('algorithms/rpca/IALM');
    warning('on','all');
  end
  %
  % APG (Lin et al. 2009)
  % process_video('RPCA', 'APG', 'dataset/demo.avi', 'output/demo_APG.avi');
  if(strcmp(algorithm_id,'APG'))
    addpath('algorithms/rpca/APG');
    timerVal = tic;
    lambda = 1/sqrt(max(size(M,1),size(M,2)));
    [L, S] = proximal_gradient_rpca(M,lambda); % imagesc(L); imagesc(S);
    cputime = toc(timerVal);
    clear lambda;
    rmpath('algorithms/rpca/APG');
  end
  %
  % APG Partial (Lin et al. 2009)
  % process_video('RPCA', 'APG_PARTIAL', 'dataset/demo.avi', 'output/demo_APG_PARTIAL.avi');
  if(strcmp(algorithm_id,'APG_PARTIAL'))
    addpath('algorithms/rpca/APG');
    addpath('libs/PROPACK');
    timerVal = tic;
    lambda = 1/sqrt(max(size(M,1),size(M,2)));
    [L, S] = partial_proximal_gradient_rpca(M,lambda); % imagesc(L); imagesc(S);
    cputime = toc(timerVal);
    clear lambda;
    rmpath('libs/PROPACK');
    rmpath('algorithms/rpca/APG');
  end
  %
  % Dual RPCA (Lin et al. 2009)
  % process_video('RPCA', 'DUAL', 'dataset/demo.avi', 'output/demo_DUAL.avi');
  if(strcmp(algorithm_id,'DUAL')) 
    addpath('algorithms/rpca/Dual Method');
    timerVal = tic;
    lambda = 1/sqrt(max(size(M,1),size(M,2)));
    [L, S] = dual_rpca_2(M,lambda); % imagesc(L); imagesc(S);
    cputime = toc(timerVal);
    clear lambda;
    rmpath('algorithms/rpca/Dual Method');
  end
  %
  % SVT (Cai et al. 2008)
  % process_video('RPCA', 'SVT', 'dataset/demo.avi', 'output/demo_SVT.avi');
  if(strcmp(algorithm_id,'SVT')) 
    addpath('algorithms/rpca/SVT');
    timerVal = tic;
    lambda = 1/sqrt(max(size(M,1),size(M,2)));
    [L, S] = singular_value_rpca(M,lambda); % imagesc(L); imagesc(S);
    cputime = toc(timerVal);
    clear lambda;
    rmpath('algorithms/rpca/SVT');
  end
  %
  % ADM / LRSD (Yuan and Yang 2009)
  % process_video('RPCA', 'ADM', 'dataset/demo.avi', 'output/demo_ADM.avi');
  % TODO: ----> Works only on win32 (mexsvd.mexw32)
  if(strcmp(algorithm_id,'ADM')) 
    addpath('algorithms/rpca/ADM');
    t = 0.01;
    opts = [];
    opts.beta = .25/mean(abs(M(:))); % 0.10;
    opts.tol = 1e-6;
    opts.maxit = 1000;
    opts.print = 1;
    timerVal = tic;
    out = ADM(M, t/(1-t), opts);
    cputime = toc(timerVal);
    L = out.LowRank; % imagesc(L);
    S = out.Sparse; % imagesc(S);
    rmpath('algorithms/rpca/ADM');
  end
  %
  % LSADM (Goldfarb et al. 2010)
  % process_video('RPCA', 'LSADM', 'dataset/demo.avi', 'output/demo_LSADM.avi');
  if(strcmp(algorithm_id,'LSADM')) 
    addpath('algorithms/rpca/LSADM');
    opts.D = M;
    opts.mu = norm(M)/1.25;
    [n1,n2] = size(M);
    opts.Xs = M;
    opts.Ys = M;
    opts.n1 = n1;
    opts.n2 = n2;
    opts.sigma = 1e-6;
    opts.maxitr = 500;
    opts.rho = 1/sqrt(n1); 
    opts.eta_mu = 2/3;
    opts.eta_sigma = 2/3;
    opts.muf = 1e-6;
    opts.sigmaf = 1e-6;
    opts.epsilon = 1e-7;
    opts.sv = 100;
    timerVal = tic;
    out_ALM = ALM_SADAL_smoothed(opts.D,opts);
    cputime = toc(timerVal);
    L = out_ALM.X; % imagesc(L);
    S = out_ALM.Y; % imagesc(S);
    rmpath('algorithms/rpca/LSADM');
  end
  %
  % L1 Filtering (Liu et al. 2011)
  % process_video('RPCA', 'L1F', 'dataset/demo.avi', 'output/demo_L1F.avi');
  if(strcmp(algorithm_id,'L1F')) 
    addpath('algorithms/rpca/L1F');
    timerVal = tic;
    [L, S] = rpca_l1f(M); % imagesc(L); imagesc(S);
    cputime = toc(timerVal);
    rmpath('algorithms/rpca/L1F');
  end
  %
  % DECOLOR (Zhou et al. 2011)
  % process_video('RPCA', 'DECOLOR', 'dataset/demo.avi', 'output/demo_DECOLOR.avi');
  % TODO: ----> Don't works in Matlab 2014a 64bits
  if(strcmp(algorithm_id,'DECOLOR')) 
    addpath('algorithms/rpca/DECOLOR');
    addpath('algorithms/rpca/DECOLOR/internal');
    addpath(genpath('algorithms/rpca/DECOLOR/gco-v3.0'));
    timerVal = tic;
    [L, S] = DECOLOR(M); % imagesc(L); imagesc(S);
    cputime = toc(timerVal);
    rmpath(genpath('algorithms/rpca/DECOLOR/gco-v3.0'));
    rmpath('algorithms/rpca/DECOLOR/internal');
    rmpath('algorithms/rpca/DECOLOR');
  end
  %
  % GoDec (Zhou and Tao 2011)
  % process_video('RPCA', 'GoDec', 'dataset/cctv.avi', 'output/demo_GoDec.avi');
  if(strcmp(algorithm_id,'GoDec'))
    addpath('algorithms/rpca/GoDec');
    rank = 1;
    card = numel(M); %card = 3.1e+5;
    power = 0;
    timerVal = tic;
    [L, S] = GoDec(M,rank,card,power); % imagesc(L); imagesc(S);
    cputime = toc(timerVal);
    clear rank power card;
    rmpath('algorithms/rpca/GoDec');
  end
  %
  % Semi-Soft GoDec (Zhou and Tao 2011)
  % process_video('RPCA', 'SSGoDec', 'dataset/cctv.avi', 'output/demo_SSGoDec.avi');
  if(strcmp(algorithm_id,'SSGoDec'))
    addpath('algorithms/rpca/GoDec');
    rank = 1;
    tau = 8;
    power = 0;
    timerVal = tic;
    L = SSGoDec(M,rank,tau,power); % imagesc(L);
    cputime = toc(timerVal);
    S = M - L; % imagesc(S);
    clear tau power rank;
    rmpath('algorithms/rpca/GoDec');
  end
  %
  % NSA v1 (Aybat et al. 2011)
  % process_video('RPCA', 'NSA1', 'dataset/demo.avi', 'output/demo_NSA1.avi');
  if(strcmp(algorithm_id,'NSA1'))
    addpath('algorithms/rpca/NSA/NSA_v1');
    addpath('algorithms/rpca/NSA/NSA_v1/Subroutines');
    addpath('libs/PROPACK_SVT');
    stdev = 1;
    tol = 5e-6; % optimality tolerance for stopping_type 1 
    timerVal = tic;
    L = nsa_v1(M,stdev,tol,1); % imagesc(L);
    cputime = toc(timerVal);
    S = M - L; % imagesc(S);
    clear stdev tol;
    rmpath('libs/PROPACK_SVT');
    rmpath('algorithms/rpca/NSA/NSA_v1/Subroutines');
    rmpath('algorithms/rpca/NSA/NSA_v1');
  end
  %
  % NSA v2 (Aybat et al. 2011)
  % process_video('RPCA', 'NSA2', 'dataset/demo.avi', 'output/demo_NSA2.avi');
  if(strcmp(algorithm_id,'NSA2'))
    addpath('algorithms/rpca/NSA/NSA_v2');
    addpath('libs/PROPACK_SVT');
    stdev = 1;
    tol = 5e-6; % optimality tolerance for stopping_type 1 
    timerVal = tic;
    [L, S] = nsa_v2(M,stdev,tol,1); % imagesc(L); imagesc(S);
    cputime = toc(timerVal);
    clear stdev tol;
    rmpath('libs/PROPACK_SVT');
    rmpath('algorithms/rpca/NSA/NSA_v2');
  end
  %
  % PSPG (Aybat et al. 2012)
  % process_video('RPCA', 'PSPG', 'dataset/demo.avi', 'output/demo_PSPG.avi');
  if(strcmp(algorithm_id,'PSPG'))
    addpath('algorithms/rpca/PSPG');
    addpath('algorithms/rpca/PSPG/Subroutines');
    addpath('libs/PROPACK_SVT');
    stdev = 1;
    tol = 0.05;
    timerVal = tic;
    L = pspg(M,stdev,tol); % imagesc(L);
    cputime = toc(timerVal);
    S = M - L; % imagesc(S);
    clear stdev tol;
    rmpath('libs/PROPACK_SVT');
    rmpath('algorithms/rpca/PSPG/Subroutines');
    rmpath('algorithms/rpca/PSPG');
  end
  %
  % Bayesian Robust PCA with Markov Dependency (Ding et al. 2011)
  % process_video('RPCA', 'BRPCA-MD', 'dataset/demo.avi', 'output/demo_BRPCA-MD.avi');
  %
  % BRPCA-MD with Non-Stationary Noise (Ding et al. 2011)
  % process_video('RPCA', 'BRPCA-MD-NSS', 'dataset/demo.avi', 'output/demo_BRPCA-MD-NSS.avi');
  if(strcmp(algorithm_id,'BRPCA-MD') || strcmp(algorithm_id,'BRPCA-MD-NSS'))
    K = 20;
    if(strcmp(algorithm_id,'BRPCA-MD'))
      addpath('algorithms/rpca/BRPCA/MarkovDependency');
      Theta0 = InitialPara_random_MarkovDep(M,K);
    end
    if(strcmp(algorithm_id,'BRPCA-MD-NSS'))
      addpath('algorithms/rpca/BRPCA/NonstationaryNoise');
      Theta0 = InitialPara_random_MarkovDep_NN(M,K);
    end
    [~,N] = size(M);
    hyperpara.a0 = 1/K;
    hyperpara.b0 = 1-hyperpara.a0;
    hyperpara.c0 = 1e-6;
    hyperpara.d0 = 1e-6;
    hyperpara.e0 = 1e-6;
    hyperpara.f0 = 1e-6;
    hyperpara.g0 = 1e-6;
    hyperpara.h0 = 1e-6;    
    hyperpara.alpha0 = 0.01*N;
    hyperpara.beta0 = 0.99*N;
    hyperpara.alpha1 = 0.99*N;
    hyperpara.beta1 = 0.01*N; 
    MCMCpara.nBurnin = 500;
    MCMCpara.nCollect = 100;
    %MCMCpara.nBurnin = 50;
    %MCMCpara.nCollect = 10;
    timerVal = tic;
    if(strcmp(algorithm_id,'BRPCA-MD'))
      out = Bayesian_RPCAmcmc_MarkovDep(M,Theta0,opts.rows,opts.cols,hyperpara,MCMCpara);
      rmpath('algorithms/rpca/BRPCA/MarkovDependency');
    end
    if(strcmp(algorithm_id,'BRPCA-MD-NSS'))
      out = Bayesian_RPCAmcmc_MarkovDep_NN(M,Theta0,opts.rows,opts.cols,hyperpara,MCMCpara);
      rmpath('algorithms/rpca/BRPCA/NonstationaryNoise');
    end
    cputime = toc(timerVal);
    L = out.Lowrank_mean; % imagesc(L);
    S = out.Sparse_mean; % imagesc(S);
    clear hyperpara MCMCpara Theta0 K out N;
  end
  %
  % Variational Bayesian RPCA (Babacan et al. 2011)
  % process_video('RPCA', 'VBRPCA', 'dataset/demo.avi', 'output/demo_VBRPCA.avi');
  if(strcmp(algorithm_id,'VBRPCA'))
    addpath('algorithms/rpca/VBRPCA');
    % all options are *optional*, everything will be set automatically
    % you can modify these options to get better performance
    options.verbose = 1;
    options.initial_rank = 'auto'; % This sets to the maximum possible rank
    % options.initial_rank = 300; % or we can use a value. 
    %options.X_true = X_true;
    %options.E_true = E_true;
    options.inf_flag = 2; % inference flag for the sparse component
    % 1 for standard VB, 2 for MacKay. MacKay generally converges faster.
    options.MAXITER = 200;
    %Estimate noise variance? (beta is inverse noise variance)
    options.UPDATE_BETA = 1; 
    % If the noise inv. variance is not to be estimated, set
    % options.UPDATE_BETA = 0; % and set beta using
    % options.beta = 1e3; 
    % Select the optimization mode: 
    % 'VB': fully Bayesian inference (default)
    % 'VB_app': fully Bayesian with covariance approximation
    % 'MAP': maximum a posteriori (covariance is set to 0)
    options.mode = 'VB';
    % For large scale problems, set to 'VB_app'. 
    % options.mode = 'VB_app';
    timerVal = tic;
    [L, ~, ~, S] = VBRPCA(M, options); % imagesc(L); % imagesc(S);
    cputime = toc(timerVal);
    clear options;
    rmpath('algorithms/rpca/VBRPCA');
  end
  %
  % PRMF (Wang et al. 2012)
  % process_video('RPCA', 'PRMF', 'dataset/demo.avi', 'output/demo_PRMF.avi');
  if(strcmp(algorithm_id,'PRMF'))
    addpath('algorithms/rpca/PRMF');
    X = normalize(M);
    rk = 2;
    lambdaU = 1;
    lambdaV = 1;
    tol = 1e-2;
    timerVal = tic;
    [P, Q] = RPMF(X, rk, lambdaU, lambdaV, tol);
    cputime = toc(timerVal);
    L = P * Q; % imagesc(L);
    %S = abs(X - P * Q);
    S = X - L; % imagesc(S);
    clear X P Q rk lambdaU lambdaV tol;
    rmpath('algorithms/rpca/PRMF');
  end
  %
  % Online PRMF (Wang et al. 2012)
  % process_video('RPCA', 'OPRMF', 'dataset/demo.avi', 'output/demo_OPRMF.avi');
  if(strcmp(algorithm_id,'OPRMF'))
    addpath('algorithms/rpca/PRMF');
    X = normalize(M);
    rk = 2;
    lambdaU = 1;
    lambdaV = 1;
    tol = 1e-2;
    mask = ones(size(X));
    timerVal = tic;
    [~, ~, L] = onlineRPMF(X, rk, lambdaU, lambdaV, tol, mask); % imagesc(L);
    cputime = toc(timerVal);
    S = X - L; % imagesc(S);
    clear X rk lambdaU lambdaV tol mask;
    rmpath('algorithms/rpca/PRMF');
  end
  %
  % Markov BRMF (Wang and Yeung 2013)
  % process_video('RPCA', 'MBRMF', 'dataset/demo.avi', 'output/demo_MBRMF.avi');
  if(strcmp(algorithm_id,'MBRMF'))
    addpath('algorithms/rpca/BRMF');
    addpath('algorithms/rpca/BRMF/Utilities');
    addpath('algorithms/rpca/BRMF/mex');
    D = normalize(M);
    % Set up the (hyper)parameters. See more details in the paper.
    % r = 10;
    r = round(min(20,sqrt(size(D,2)))/2);
    opts.maxIter = 100;
    opts.burnin = 50;
    opts.invW_0 = 1000 * eye(r * 2);
    opts.beta_0 = 2;
    opts.nu_0 = r * 2;
    opts.a = 1e-4;
    opts.b = 1e0; %[1 ~ 10]
    % We set the maximum rank to be twice of the ground truth.
    opts.r = r * 2;
    opts.alpha = 0.5;
    timerVal = tic;
    [~, ~, ~, p] = MBRMF(D, opts);
    cputime = toc(timerVal);
    L = p; % imagesc(L);
    S = (D - p); % imagesc(S);
    clear D r opts p;
    rmpath('algorithms/rpca/BRMF');
    rmpath('algorithms/rpca/BRMF/Utilities');
    rmpath('algorithms/rpca/BRMF/mex');
  end
  %
  % TFOCS (Becker et al. 2011)
  % process_video('RPCA', 'TFOCS-EC', 'dataset/demo.avi', 'output/demo_TFOCS-EC.avi');
  % process_video('RPCA', 'TFOCS-IC', 'dataset/demo.avi', 'output/demo_TFOCS-IC.avi');
  % video = load_video_file('dataset/cctv.avi');
  if(strcmp(algorithm_id,'TFOCS-EC') || strcmp(algorithm_id,'TFOCS-IC'))
    addpath('algorithms/rpca/TFOCS');
    addpath('libs/TFOCS-1.3.1');
    timerVal = tic;
    if(strcmp(algorithm_id,'TFOCS-EC'))
      L = tfocs_interface(M, 1); % imagesc(L); imagesc(S);
    end
    if(strcmp(algorithm_id,'TFOCS-IC'))
      L = tfocs_interface(M, 2);
    end
    cputime = toc(timerVal);
    S = M - L;
    rmpath('libs/TFOCS-1.3.1');
    rmpath('algorithms/rpca/TFOCS');
  end
  %
  % A variational approach to SPCP (Aravkin et al. 2014)
  %
  % RPCA | SPCP-sum-SPG | Stable PCP-sum solved by Spectral Projected Gradient (Aravkin et al. 2014)
  % RPCA | SPCP-max-QN  | Stable PCP-max solved by Quasi-Newton (Aravkin et al. 2014)
  % RPCA | Lag-SPCP-SPG | Lagrangian SPCP solved by Spectral Projected Gradient (Aravkin et al. 2014)
  % RPCA | Lag-SPCP-QN  | Lagrangian SPCP solved by Quasi-Newton (Aravkin et al. 2014)
  %
  % process_video('RPCA', 'flip-SPCP-sum-SPG', 'dataset/demo.avi', 'output/demo_flip-SPCP-sum-SPG.avi');
  % process_video('RPCA', 'flip-SPCP-max-QN', 'dataset/demo.avi', 'output/demo_flip-SPCP-max-QN.avi');
  % process_video('RPCA', 'Lag-SPCP-SPG', 'dataset/demo.avi', 'output/demo_Lag-SPCP-SPG.avi');
  % process_video('RPCA', 'Lag-SPCP-QN', 'dataset/demo.avi', 'output/demo_Lag-SPCP-QN.avi');
  if(strcmp(algorithm_id,'flip-SPCP-sum-SPG') ...
  || strcmp(algorithm_id,'flip-SPCP-max-QN') ...
  || strcmp(algorithm_id,'Lag-SPCP-SPG') ...
  || strcmp(algorithm_id,'Lag-SPCP-QN')) 
    addpath('algorithms/rpca/SPGL1');
    nFrames     = size(M,2);
    lambda      = 1/sqrt(max(size(M,1),size(M,2)));
    L0          = repmat(median(M,2), 1, nFrames);
    S0          = M - L0;
    epsilon     = 5e-3*norm(M,'fro'); % tolerance for fidelity to data
    
    timerVal = tic;
    if(strcmp(algorithm_id,'flip-SPCP-sum-SPG')) % Flip-Flop version of SPCP-sum solved by Spectral Projected Gradient
      opts = struct('sum',true,'L0',L0,'S0',S0,'max',false,...
                    'tau0',3e5,'SPGL1_tol',1e-1,'tol',1e-3);
      [L,S] = solver_RPCA_SPGL1(M,lambda,epsilon,[],opts);
    end
    if(strcmp(algorithm_id,'flip-SPCP-max-QN')) % Flip-Flop version pf SPCP-max solved by Quasi-Newton
      opts = struct('sum',false,'L0',L0,'S0',S0,'max',true,...
                    'tau0',3e5,'SPGL1_tol',1e-1,'tol',1e-3);
      [L,S] = solver_RPCA_SPGL1(M,lambda,epsilon,[],opts);
    end
    if(strcmp(algorithm_id,'Lag-SPCP-SPG')) % Lagrangian SPCP solved by Spectral Projected Gradient
      opts    = struct('sum',false,'L0',L0,'S0',S0,'max',false,'tol',1e-3);
      lambdaL = 0.25; lambdaS = 0.01; % (Aravkin et al. 2014)
      [L,S] = solver_RPCA_Lagrangian(M,lambdaL,lambdaS,[],opts);
    end
    if(strcmp(algorithm_id,'Lag-SPCP-QN')) % Lagrangian SPCP solved by Quasi-Newton
      opts    = struct('sum',false,'L0',L0,'S0',S0,'max',false,'tol',1e-3,'quasiNewton',true);
      lambdaL = 0.25; lambdaS = 0.01;
      [L,S] = solver_RPCA_Lagrangian(M,lambdaL,lambdaS,[],opts);
    end
    cputime = toc(timerVal);
    
    clear nFrames lambda lambdaL lambdaS L0 S0 epsilon opts;
    rmpath('algorithms/rpca/SPGL1');
  end
  %
  %
  results.L = L; % low-rank matrix
  results.S = S; % sparse matrix
  results.cputime = cputime;
  clear L S;
end
