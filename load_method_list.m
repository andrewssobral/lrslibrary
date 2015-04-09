%% cell load_method_list(void)
%
function [L] = load_method_list()
  fileID = fopen('lrs_methods.dat');
  L = textscan(fileID,'%s %s', 'delimiter', '|');
  fclose(fileID);
end
