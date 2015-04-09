%% [struct] = process_tensor(string, string, tensor);
%
function [results] = process_tensor(method_id, algorithm_id, T)
  method_name = get_method_name_by_id(method_id);
  algorithm_name = get_algorithm_name_by_id(algorithm_id);
  displog(['Running ' method_name ' with ' algorithm_name]);

  %%% NTF methods
  % i.e: results = process_tensor('NTF', 'betaNTF', T);
  if(strcmp(method_id,'NTF'))
    results = run_algorithm_ntf(algorithm_id, T);
  end
  
  %%% TD methods
  % i.e: results = process_tensor('TD', 'HOSVD', T);
  if(strcmp(method_id,'TD'))
    results = run_algorithm_td(algorithm_id, T);
  end
  
  %%% Apply hard thresholding
  results.O = hard_threshold(results.S);
  
  displog('Decomposition finished');
end
