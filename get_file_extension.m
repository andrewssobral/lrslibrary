%% string get_file_extension(string)
% filePath - string
% ext - string
% 
function [ext] = get_file_extension(fileName)
  ext = fileName(length(fileName)-2:length(fileName));
end

