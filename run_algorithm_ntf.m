%%% NTF algorithms
% struct = run_algorithm_ntf(string, tensor)
%
function results = run_algorithm_ntf(algorithm_id, T)
  load_lrslibrary_conf;
  
  alg_path = fullfile(lrs_conf.ntf_path,algorithm_id);
  addpath(genpath(alg_path));
  
  L = []; % low-rank tensor
  S = []; % sparse tensor
  results.cputime = 0;
  
  timerVal = tic;
  % warning('off','all');
  
  %
  % NTF | betaNTF | Simple beta-NTF implementation (Antoine Liutkus, 2012)
  % process_video('NTF', 'betaNTF', 'dataset/demo.avi', 'output/demo_beta-NTF.avi');
  %
  if(strcmp(algorithm_id,'betaNTF'))
    % Compute a simple NTF model of 10 components
    A = double(T);
    r = 10;
    [~,~,~,L] = betaNTF(A,r);
    S = (A - L);
    % For reconstruction
    % for i = 1:r, B_hat(:,:,i) = W * diag(Q(i,:)) * H'; end
  end
  %
  % NTF | NTD-MU   | Non-negative Tucker Decomposition solved by Multiplicative Updates (Zhou et al. 2012)
  % NTF | NTD-APG  | Non-negative Tucker Decomposition solved by Accelerated Proximal Gradient (Zhou et al. 2012)
  % NTF | NTD-HALS | Non-negative Tucker Decomposition solved by Hierarchical ALS  (Zhou et al. 2012)
  %
  % process_video('NTF', 'NTD-MU', 'dataset/demo.avi', 'output/demo_NTD-MU.avi');
  % process_video('NTF', 'NTD-APG', 'dataset/demo.avi', 'output/demo_NTD-APG.avi');
  % process_video('NTF', 'NTD-HALS', 'dataset/demo.avi', 'output/demo_NTD-HALS.avi');
  %
  if(strcmp(algorithm_id,'NTD-MU') ...
  || strcmp(algorithm_id,'NTD-APG') ...
  || strcmp(algorithm_id,'NTD-HALS'))
    alg_path = fullfile(lrs_conf.ntf_path,'lraNTD');
    addpath(genpath(alg_path));
    
    R = [size(T,1) size(T,2) 2];
    if(strcmp(algorithm_id,'NTD-MU'))
      opts = struct('NumOfComp',R,'nlssolver','mu','maxiter',100,'maxiniter',20,'tdalgFile','call_tucker_als_opts.mat');
    end
    if(strcmp(algorithm_id,'NTD-APG'))
      opts = struct('NumOfComp',R,'nlssolver','apg','maxiter',100,'maxiniter',20,'tdalgFile','call_tucker_als_opts.mat');
    end
    if(strcmp(algorithm_id,'NTD-HALS'))
      opts = struct('NumOfComp',R,'nlssolver','hals','maxiter',100,'maxiniter',20,'tdalgFile','call_tucker_als_opts.mat');
    end
    
    T_hat = lraNTD_ANLS(T, opts);
    
    L = double(tensor(T_hat));
    S = double(T) - L;
  end
  %
  % NTF | bcuNTD | Non-negative Tucker Decomposition by block-coordinate update (Xu and Yin, 2012)
  % process_video('NTF', 'bcuNTD', 'dataset/demo.avi', 'output/demo_bcuNTD.avi');
  %
  if(strcmp(algorithm_id,'bcuNTD'))
    % Compute a simple NTF model of 10 components
    A = double(T);
    R = [size(T,1) size(T,2) 2]; % tensor rank
    opts.maxit = 1000; % max number of iterations
    opts.tol = 1e-4; % stopping tolerance
    [M,C] = ntd(T,R,opts);
    L = (double(full(ttensor(C,M))));
    S = (A - L);
  end
  %
  % NTF | bcuNCP | Non-negative CP Decomposition by block-coordinate update (Xu and Yin, 2012)
  % process_video('NTF', 'bcuNCP', 'dataset/demo.avi', 'output/bcuNCP.avi');
  %
  if(strcmp(algorithm_id,'bcuNCP'))
    % Compute a simple NTF model of 10 components
    A = double(T);
    R = 10; % tensor rank
    opts.maxit = 1000; % max number of iterations
    opts.tol = 1e-4; % stopping tolerance
    M = ncp(T,R,opts);
    L = double(full(M));
    S = (A - L);
  end
  %
  cputime = toc(timerVal);
  rmpath(genpath(alg_path));
  %
  results.L = L; % low-rank tensor
  results.S = S; % sparse tensor
  results.cputime = cputime;
end
