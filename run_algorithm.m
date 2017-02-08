%%% Launch algorithms
% results = run_algorithm(method_id, algorithm_id, A, params)
%
% method_id - string
% algorithm_id - string
% A - matrix/tensor
% params - struct (optional)
%
% params.Idx - observed indexes (vector)
% params.Omega - observed entries (binary matrix/tensor)
%
% results - struct
% results.L - low-rank component
% results.S - sparse component
% results.O - hard thresholding
% results.cputime - CPU time;
%
function results = run_algorithm(method_id, algorithm_id, data, params)
  if(nargin < 3)
    error('method, algorithm and data must be defined');
  end
  if(nargin < 4)
    params = [];
  end
  lrs_load_conf;
  
  method_id = strtrim(method_id);
  algorithm_id = strtrim(algorithm_id);
  
  method_name = get_method_name_by_id(method_id);
  algorithm_name = get_algorithm_name_by_id(algorithm_id);
  displog(['Running ' method_name ' with ' algorithm_name]);
	
  switch method_id
		case 'RPCA'
			method_path = lrs_conf.rpca_path;
		case 'ST'
			method_path = lrs_conf.st_path;
		case 'MC'
			method_path = lrs_conf.mc_path;
		case 'LRR'
			method_path = lrs_conf.lrr_path;
		case 'TTD'
			method_path = lrs_conf.ttd_path;
		case 'NMF'
			method_path = lrs_conf.nmf_path;
		case 'NTF'
			method_path = lrs_conf.ntf_path;
		case 'TD'
			method_path = lrs_conf.td_path;
		otherwise
			error('Undefined method!');
	end
  
  alg_path = fullfile(method_path,algorithm_id);
  addpath(genpath(alg_path));
  
  params.A = data; M = data; T = data;
  L = zeros(size(data)); % low-rank component
  S = zeros(size(data)); % sparse component
  results.cputime = 0;
  
  if(~isfield(params,'rows') && ~isfield(params,'cols'))
    params.rows = size(data,1);
    params.cols = size(data,2);
  end
  
  % For matrix/tensor completion
  if(~isfield(params,'Idx') && ~isfield(params,'Omega')) 
    [params.Idx, params.Omega] = subsampling(data, 0.5);
  end
  if(isfield(params,'Idx') && ~isfield(params,'Omega')) % Build "params.Omega" from "params.Idx"
    params.Omega = zeros(size(data)); 
    params.Omega(params.Idx) = 1;
  end
  if(isfield(params,'Omega')) % Build "params.Idx" from "params.Omega"
    % params.Omega = ones(size(A)); 
    % params.Omega = randi([0 1],size(A)); 
    params.Idx = find(params.Omega);
  end
  Idx = params.Idx; % MIdx = M(Idx);
  Omega = params.Omega; % MOmega = M.*Omega;
  
  timerVal = tic;
  % warning('off','all');
  run_alg;
  % warning('on','all');
  cputime = toc(timerVal);
  
  rmpath(genpath(alg_path));
  
  results.L = L; % low-rank component
  results.S = S; % sparse component
  results.O = hard_threshold(S); % apply hard thresholding
  results.Omega = params.Omega; % for matrix/tensor completion
  results.Idx = params.Idx; % for matrix/tensor completion
  results.cputime = cputime;
	displog('Decomposition finished');
end
