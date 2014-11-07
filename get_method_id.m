%% int get_method_id(string)
% method - string
% id - int
function [method_id] = get_method_id(method)
  if(length(method) > 1)
    method_split = strsplit(method,'-');
    method_id = strtrim(method_split(1));
    method_id = char(method_id);
  else
    method_id = 0;
  end
end
