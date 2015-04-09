%% string get_method_name_by_id(int)
% method_id - int
% method_name - string
%
function [method_name] = get_method_name_by_id(method_id)
  %clear all; clc; method_id = 'RPCA';
  method_name = '-';
  L = load_method_list();
  [m,n] = size(L{1,1});
  for i = 1 : m
    m_id = strtrim(char(L{1,1}(i)));
    m_name = strtrim(char(L{1,2}(i)));
    
    if(strcmp(method_id,m_id))
      method_name = m_name;
      break;
    end
  end
end

