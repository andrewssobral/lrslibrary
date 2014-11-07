%% cell load_algorithm_list(void)
%
function [L] = load_algorithm_list()
  fileID = fopen('algorithms.dat');
  L = textscan(fileID,'%s %s %s %d', 'delimiter', '|');
  fclose(fileID);
end
