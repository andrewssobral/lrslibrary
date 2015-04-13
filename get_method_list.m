%% cell fill_method_listbox(void)
%
function [list] = get_method_list()
  L = load_method_list();
  [m,n] = size(L{1,1});
  list = {};
  k = 1;
  for i = 1 : m
    method_id = strtrim(char(L{1,1}(i)));
    method_name = strtrim(char(L{1,2}(i)));
    
    if(method_id == '-')
      method = [method_id];
    else
      method = [method_id ' - ' method_name];
    end
    
    list{k,1} = method;
    k = k + 1;
  end
end

