%% string get_algorithm_name_by_number(int)
% algorithm_id - int
% algorithm_name - string
%
function [algorithm_name] = get_algorithm_name_by_id(algorithm_id)
  %clear all; clc; algorithm_id = 'FPCP';
  L = load_algorithm_list();
  [m,n] = size(L{1,1});
  algorithm_name = '';
  for i = 1 : m
    a_id = strtrim(char(L{1,2}(i)));
    a_name = strtrim(char(L{1,3}(i)));
    %disp([algorithm_id '<->' a_id]);
    if(strcmp(algorithm_id,a_id))
      algorithm_name = a_name;
      break;
    end
  end
end

