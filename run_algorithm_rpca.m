%%% RPCA algorithms
% struct = run_algorithm_rpca(string, 2dmatrix)
%
function results = run_algorithm_rpca(algorithm_id, M, opts)
  lrs_load_conf;
  
  alg_path = fullfile(lrs_conf.rpca_path,algorithm_id);
  addpath(genpath(alg_path));
  
  lambda = 1/sqrt(max(size(M))); % default lambda
  
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
    % RPCA (De la Torre and Black, 2001)
    % process_video('RPCA', 'RPCA', 'dataset/demo.avi', 'output/demo_RPCA.avi');
    if(strcmp(algorithm_id,'RPCA'))
      sizeim = [opts.rows opts.cols];
      [L,S] = run_RPCA(M,sizeim);
    end
    %
    % PCP (Candes et al. 2009)
    % process_video('RPCA', 'PCP', 'dataset/demo.avi', 'output/demo_PCP.avi');
    if(strcmp(algorithm_id,'PCP'))
      tol = 1e-5;
      [L,S] = PCP(M,lambda,tol);
    end
    %
    % Fast PCP (Rodriguez and Wohlberg 2013)
    % process_video('RPCA', 'FPCP', 'dataset/demo.avi', 'output/demo_FPCP.avi');
    if(strcmp(algorithm_id,'FPCP'))
      [L,S] = fastpcp(M,lambda);
    end
    %
    % R2PCP: Riemannian Robust Principal Component Pursuit (Hintermüller and Wu, 2014)
    % process_video('RPCA', 'R2PCP', 'dataset/demo.avi', 'output/demo_R2PCP.avi');
    if(strcmp(algorithm_id,'R2PCP'))
      run_R2PCP;
    end
    %
    % AS-RPCA: Active Subspace: Towards Scalable Low-Rank Learning (Liu and Yan, 2012)
    % process_video('RPCA', 'AS-RPCA', 'dataset/demo.avi', 'output/demo_AS-RPCA.avi');
    if(strcmp(algorithm_id,'AS-RPCA'))
      lambda = 1/sqrt(min(size(M)));
      [L,S] = as_rpca(M,lambda);
    end
    %
    % ALM (Tang and Nehorai 2011)
    % process_video('RPCA', 'ALM', 'dataset/demo.avi', 'output/demo_ALM.avi');
    if(strcmp(algorithm_id,'ALM'))
      [L,S] = alm(M);
    end
    %
    % EALM (Lin et al. 2009)
    % process_video('RPCA', 'EALM', 'dataset/demo.avi', 'output/demo_EALM.avi');
    if(strcmp(algorithm_id,'EALM'))
      [L,S] = exact_alm_rpca(M);
    end
    %
    % IALM (Lin et al. 2009)
    % process_video('RPCA', 'IALM', 'dataset/demo.avi', 'output/demo_IALM.avi');
    if(strcmp(algorithm_id,'IALM'))
      [L,S] = inexact_alm_rpca(M);
    end
    %
    % IALM + LMSVDS (Liu et al. 2012)
    % process_video('RPCA', 'IALM_LMSVDS', 'dataset/demo.avi', 'output/demo_IALM_LMSVDS.avi');
    if(strcmp(algorithm_id,'IALM_LMSVDS'))
      [L,S] = inexact_alm_rpca_with_lmsvds(M);
    end
    %
    % IALM + BLWS (Lin and Wei 2010)
    % process_video('RPCA', 'IALM_BLWS', 'dataset/demo.avi', 'output/demo_IALM_BLWS.avi');
    if(strcmp(algorithm_id,'IALM_BLWS'))
      [L,S] = inexact_alm_rpca_with_blws(M);
    end
    %
    % APG (Lin et al. 2009)
    % process_video('RPCA', 'APG', 'dataset/demo.avi', 'output/demo_APG.avi');
    if(strcmp(algorithm_id,'APG'))
      [L,S] = proximal_gradient_rpca(M,lambda);
    end
    %
    % APG Partial (Lin et al. 2009)
    % process_video('RPCA', 'APG_PARTIAL', 'dataset/demo.avi', 'output/demo_APG_PARTIAL.avi');
    if(strcmp(algorithm_id,'APG_PARTIAL'))
      [L,S] = partial_proximal_gradient_rpca(M,lambda);
    end
    %
    % Dual RPCA (Lin et al. 2009)
    % process_video('RPCA', 'DUAL', 'dataset/demo.avi', 'output/demo_DUAL.avi');
    if(strcmp(algorithm_id,'DUAL'))
      [L,S] = dual_rpca_2(M,lambda);
    end
    %
    % SVT (Cai et al. 2008)
    % process_video('RPCA', 'SVT', 'dataset/demo.avi', 'output/demo_SVT.avi');
    if(strcmp(algorithm_id,'SVT'))
      [L,S] = singular_value_rpca(M,lambda,1e4,0.9,'svd');
    end
    %
    % ADM / LRSD (Yuan and Yang 2009)
    % process_video('RPCA', 'ADM', 'dataset/demo.avi', 'output/demo_ADM.avi');
    % TODO: ----> Works only on win32 (mexsvd.mexw32)
    if(strcmp(algorithm_id,'ADM')) 
      run_ADM;
    end
    %
    % LSADM (Goldfarb et al. 2010)
    % process_video('RPCA', 'LSADM', 'dataset/demo.avi', 'output/demo_LSADM.avi');
    if(strcmp(algorithm_id,'LSADM')) 
      run_LSADM;
    end
    %
    % L1 Filtering (Liu et al. 2011)
    % process_video('RPCA', 'L1F', 'dataset/demo.avi', 'output/demo_L1F.avi');
    if(strcmp(algorithm_id,'L1F'))
      [L,S] = rpca_l1f(M);
    end
    %
    % DECOLOR (Zhou et al. 2011)
    % process_video('RPCA', 'DECOLOR', 'dataset/demo.avi', 'output/demo_DECOLOR.avi');
    if(strcmp(algorithm_id,'DECOLOR'))
      [L,S] = DECOLOR(M);
    end
    %
    % GoDec (Zhou and Tao 2011)
    % process_video('RPCA', 'GoDec', 'dataset/cctv.avi', 'output/demo_GoDec.avi');
    if(strcmp(algorithm_id,'GoDec'))
      rank = 1;
      card = numel(M); %card = 3.1e+5;
      power = 0;
      [L,S] = GoDec(M,rank,card,power);
    end
    %
    % Semi-Soft GoDec (Zhou and Tao 2011)
    % process_video('RPCA', 'SSGoDec', 'dataset/cctv.avi', 'output/demo_SSGoDec.avi');
    if(strcmp(algorithm_id,'SSGoDec'))
      rank = 1;
      tau = 8;
      power = 0;
      L = SSGoDec(M,rank,tau,power);
      S = M - L;
    end
    %
    % GreGoDec (Zhou and Tao 2013)
    % process_video('RPCA', 'GreGoDec', 'dataset/cctv.avi', 'output/demo_GreGoDec.avi');
    if(strcmp(algorithm_id,'GreGoDec'))
      rank = 1;
      tau = 7;
      power = 5;
      tol = 1e-3;
      k = 1;
      L = GreGoDec(M,rank,tau,tol,power,k);
      S = M - L;
    end
    %
    % NSA v1 (Aybat et al. 2011)
    % process_video('RPCA', 'NSA1', 'dataset/demo.avi', 'output/demo_NSA1.avi');
    if(strcmp(algorithm_id,'NSA1'))
      stdev = 1;
      tol = 5e-6; % optimality tolerance for stopping_type 1 
      L = nsa_v1(M,stdev,tol,1);
      S = M - L;
    end
    %
    % NSA v2 (Aybat et al. 2011)
    % process_video('RPCA', 'NSA2', 'dataset/demo.avi', 'output/demo_NSA2.avi');
    if(strcmp(algorithm_id,'NSA2'))
      stdev = 1;
      tol = 5e-6; % optimality tolerance for stopping_type 1 
      [L,S] = nsa_v2(M,stdev,tol,1);
    end
    %
    % PSPG (Aybat et al. 2012)
    % process_video('RPCA', 'PSPG', 'dataset/demo.avi', 'output/demo_PSPG.avi');
    if(strcmp(algorithm_id,'PSPG'))
      stdev = 1;
      tol = 0.05;
      L = pspg(M,stdev,tol);
      S = M - L;
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
        Theta0 = InitialPara_random_MarkovDep(M,K);
      end
      if(strcmp(algorithm_id,'BRPCA-MD-NSS'))
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
      end
      if(strcmp(algorithm_id,'BRPCA-MD-NSS'))
        out = Bayesian_RPCAmcmc_MarkovDep_NN(M,Theta0,opts.rows,opts.cols,hyperpara,MCMCpara);
      end
      L = out.Lowrank_mean;
      S = out.Sparse_mean;
    end
    %
    % Variational Bayesian RPCA (Babacan et al. 2011)
    % process_video('RPCA', 'VBRPCA', 'dataset/demo.avi', 'output/demo_VBRPCA.avi');
    if(strcmp(algorithm_id,'VBRPCA'))
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
      [L,~,~,S] = VBRPCA(M, options);
    end
    %
    % PRMF (Wang et al. 2012)
    % process_video('RPCA', 'PRMF', 'dataset/demo.avi', 'output/demo_PRMF.avi');
    if(strcmp(algorithm_id,'PRMF'))
      X = normalize(M);
      rk = 2;
      lambdaU = 1;
      lambdaV = 1;
      tol = 1e-2;
      [P, Q] = RPMF(X, rk, lambdaU, lambdaV, tol);
      L = P * Q;
      %S = abs(X - P * Q);
      S = X - L;
    end
    %
    % Online PRMF (Wang et al. 2012)
    % process_video('RPCA', 'OPRMF', 'dataset/demo.avi', 'output/demo_OPRMF.avi');
    if(strcmp(algorithm_id,'OPRMF'))
      X = normalize(M);
      rk = 2;
      lambdaU = 1;
      lambdaV = 1;
      tol = 1e-2;
      mask = ones(size(X));
      [~, ~, L] = onlineRPMF(X, rk, lambdaU, lambdaV, tol, mask);
      S = X - L;
    end
    %
    % Markov BRMF (Wang and Yeung 2013)
    % process_video('RPCA', 'MBRMF', 'dataset/demo.avi', 'output/demo_MBRMF.avi');
    if(strcmp(algorithm_id,'MBRMF'))
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
      [~,~,~,p] = MBRMF(D, opts);
      L = p;
      S = (D - p);
    end
    %
    % TFOCS (Becker et al. 2011)
    % process_video('RPCA', 'TFOCS-EC', 'dataset/demo.avi', 'output/demo_TFOCS-EC.avi');
    % process_video('RPCA', 'TFOCS-IC', 'dataset/demo.avi', 'output/demo_TFOCS-IC.avi');
    % video = load_video_file('dataset/cctv.avi');
    if(strcmp(algorithm_id,'TFOCS-EC') || strcmp(algorithm_id,'TFOCS-IC'))
      alg_path = fullfile(lrs_conf.rpca_path,'TFOCS');
      addpath(genpath(alg_path));
      if(strcmp(algorithm_id,'TFOCS-EC'))
        L = tfocs_interface(M, 1);
      end
      if(strcmp(algorithm_id,'TFOCS-IC'))
        L = tfocs_interface(M, 2);
      end
      S = M - L;
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
      alg_path = fullfile(lrs_conf.rpca_path,'SPGL1');
      addpath(genpath(alg_path));
      nFrames     = size(M,2);
      lambda      = 1/sqrt(max(size(M,1),size(M,2)));
      L0          = repmat(median(M,2), 1, nFrames);
      S0          = M - L0;
      epsilon     = 5e-3*norm(M,'fro'); % tolerance for fidelity to data
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
    end
    %
    % RPCA | RegL1-ALM | Low-Rank Matrix Approximation under Robust L1-Norm (Zheng et al. 2012)
    % process_video('RPCA', 'RegL1-ALM', 'dataset/demo.avi', 'output/demo_RegL1-ALM.avi');
    %
    if(strcmp(algorithm_id,'RegL1-ALM')) 
      W = ones(size(M));
      r = 1;
      lambda = 1e-3;
      rho = 1.2;
      maxIterIN = 1;
      signM = 0;
      %[M_est,U_est,V_est,L1_error] = ...
      L = RobustApproximation_M_UV_TraceNormReg(M,W,r,lambda,rho,maxIterIN,signM);
      S = M - L;
    end
    %
    % FW-T: SPCP solved by Frank-Wolfe method (Mu et al. 2014)
    % process_video('RPCA', 'FW-T', 'dataset/demo.avi', 'output/demo_FW-T.avi');
    if(strcmp(algorithm_id,'FW-T'))
      [m,n] = size(M); 
      D = M/norm(M,'fro'); % imagesc(M); imagesc(D);

      % parameter tuning
      rho = 1; % rho = 0.5;  % sampling ratio
      Omega = rand(m,n) <= rho; % support of observation imagesc(Omega);
      obs = Omega.*D; % measurements imagesc(obs);

      % this is parameter to control noise level
      % the smaller the noise, the smaller is delta
      delta = 0.01;

      lambda_1 = delta*rho*norm(obs,'fro'); 
      lambda_2 = delta*sqrt(rho)*norm(obs,'fro')/sqrt(max(m,n));

      par.M = D; 
      par.lambda_1 = lambda_1;
      par.lambda_2 = lambda_2;
      par.iter = 1000;
      par.display = 1;
      par.rho = rho;
      par.epsilon = 10^-3; % stopping criterion
      par.method = 'exact'; % 'exact' or 'power'
      par.Omega = Omega; % ones(m,n)
      par.compare = 0; % make comparison or not

      output = FW_T(par); % main function
      %output = fista(par); % fista function
      %output = ista(par); % ista function

      L = output.L;
      S = output.S;
    end
    %
    % GA, GM, TGA (Hauberg et al. 2014)
    % process_video('RPCA', 'TGA', 'dataset/demo.avi', 'output/demo_TGA.avi');
    if(strcmp(algorithm_id,'GA') || strcmp(algorithm_id,'GM')...
    || strcmp(algorithm_id,'TGA'))
      alg_path = fullfile(lrs_conf.rpca_path,'GA');
      addpath(genpath(alg_path));
      if(strcmp(algorithm_id,'GA'))
        L = grassmann_average(M', 1);
      end
      if(strcmp(algorithm_id,'GM'))
        L = grassmann_median(M', 1);
      end
      if(strcmp(algorithm_id,'TGA'))
        L = trimmed_grassmann_average(M', 50, 1);
      end
      L = nma_rescale(L,min(M(:)),max(M(:))); % 0.968 0.207
      L = repmat(L,1,size(M,2));
      S = M - L;
      % show_2dvideo(M,m,n);
      % show_2dvideo(L,m,n);
      % show_2dvideo(S,m,n);
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
