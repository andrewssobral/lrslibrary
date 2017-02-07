function res = getmexopts(Tag)
% function res = getmexopts(Tag)
% Get the MCC or MEX configuration
% Author Bruno Luong <brunoluong@yahoo.com>
% Last update: 29-Jun-2009

if ispc()
    optpath=prefdir;
    optfile=[optpath filesep 'compopts.bat'];
    mexoptfile=[optpath filesep 'mexopts.bat'];
else
    optpath=matlabroot;
    optfile=[optpath '/bin/mbuildopts.sh'];
    mexoptfile=[optpath '/bin/mexopts.sh']; % not sure correct path
end

% Try to get MEX option first
fid=fopen(mexoptfile,'r');
if fid<=0
    % Next MCC options
    fid=fopen(optfile,'r');
end

if fid>0
    iscompilerline=@(S) (strcmp(S,['set ' Tag]));
    C=textscan(fid,'%s %s', 'delimiter', '=', 'whitespace', '');
    fclose(fid);
    cline=find(cellfun(iscompilerline,C{1}));
    if isempty(cline)
        error('getmexopt [Bruno]: cannot get Tag %s', Tag)
    end
    res=C{2}{cline};
    root=regexprep(matlabroot,'\\','\\\\');
    res = regexprep(res,'%MATLAB%',root);
else
    error('getmexopts [Bruno]: cannot open comopts.bat file')
end

% Bruno 