%%% LRR algorithms
% struct = run_algorithm_lrr(string, 2dmatrix)
%
function results = run_algorithm(method_id, algorithm_id, data, params)
  lrs_load_conf;
  
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
  
  M = data; T = data;
  L = zeros(size(data)); % low-rank component
  S = zeros(size(data)); % sparse component
  
  results.cputime = 0;
  if(isempty(params))
    params.rows = size(data,1);
    params.cols = size(data,2);
  end
  
  timerVal = tic;
  % warning('off','all');
  run_alg;
  % warning('on','all');
  cputime = toc(timerVal);
  
  rmpath(genpath(alg_path));
  
  results.L = L; % low-rank component
  results.S = S; % sparse component
  results.O = hard_threshold(S); % apply hard thresholding
  results.cputime = cputime;
	displog('Decomposition finished');
end
