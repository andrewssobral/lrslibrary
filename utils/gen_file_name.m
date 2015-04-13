%% string gen_file_name(string,string)
% fileName - string
% extra - string
% newFileName - string
%
function [newFileName] = gen_file_name(fileName,extra)
  %ext = get_file_extension(fileName);
  %newFileName = [fileName(1:length(fileName)-4) extra ext];
  [pathstr,name,ext] = fileparts(fileName);
  newFileName = fullfile(pathstr,[name extra ext]);
end
