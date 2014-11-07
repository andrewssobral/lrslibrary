%% cell get_algorithms_list(int)
%
function [alg_list] = get_algorithm_list(id)
  L = load_algorithm_list();
  [m,~] = size(L{1,1});
  alg_list = {};
  k = 1;
  for i = 1 : m
    method_id_field = char(L{1,1}(i));
    method_id_list = strsplit(method_id_field,',');
    
    [~,s] = size(method_id_list);
    
    flag = 0;
    for j = 1 : s
      method_id = strtrim(char(method_id_list{1,j}));
      if(strcmp(method_id,'-'))
        alg_list{k,1} = '-';
        alg_list{k,2} = '-';
        alg_list{k,3} = 0;
        alg_list{k,4} = '-';
        k = k + 1;
        break;
      end
%       if(method_id(1) == '-')
%         list{i,1} = '-';
%         break;
%       end
      if(strcmp(id,method_id))
        flag = 1;
      end
    end
    
    if(flag == 0)
      continue;
    end
        
    algorithm_id = strtrim(char(L{1,2}(i)));
    algorithm_name = strtrim(char(L{1,3}(i)));
    
    algorithm = [algorithm_id ' ' algorithm_name];
    
    alg_list{k,1} = algorithm_id;
    alg_list{k,2} = algorithm_name;
    alg_list{k,3} = L{1,4}(i);
    alg_list{k,4} = algorithm;
    k = k + 1;
  end
end

