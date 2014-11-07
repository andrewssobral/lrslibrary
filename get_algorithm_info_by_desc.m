%% struct get_algorithm_info_by_desc(string,cell)
% algorithm_desc - string
% alg_list - cell
function [alg_info] = get_algorithm_info_by_desc(algorithm_desc, alg_list)
  alg_info.algorithm_id = '-';
  alg_info.algorithm_name = '-';
  alg_info.algorithm_time = 0;
  alg_info.algorithm_desc = 0;
  [r,~] = size(alg_list);
  for i = 1:r
    aux_algorithm_id = alg_list{i,1};
    aux_algorithm_name = alg_list{i,2};
    aux_algorithm_time = alg_list{i,3};
    aux_algorithm_desc = alg_list{i,4};
    if(strcmp(algorithm_desc,aux_algorithm_desc))
      alg_info.algorithm_id = aux_algorithm_id;
      alg_info.algorithm_name = aux_algorithm_name;
      alg_info.algorithm_time = aux_algorithm_time;
      alg_info.algorithm_desc = aux_algorithm_desc;
      break;
    end
  end
end

% function [algorithm_id] = get_algorithm_id(algorithm)
%   if(length(algorithm) > 1)
%     algorithm_split = strsplit(algorithm,' ');
%     algorithm_id = strtrim(algorithm_split(1));
%     algorithm_id = char(algorithm_id);
%   else
%     algorithm_id = '';
%   end
% end