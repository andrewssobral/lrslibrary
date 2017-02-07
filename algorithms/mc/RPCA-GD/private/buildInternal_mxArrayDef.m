function content = buildInternal_mxArrayDef(mxArraydefFilename)
% function content = buildInternal_mxArrayDef(mxArraydefFilename)
%
% Building the typedef of internal structure MxArray by looking inside
% the header file include file MATRIX.H. This ensure the definition used
% is compatible with the Matlab version
% The internal definition will be used by MEX file inplacecolumnmex and
% releaseinplace
%
% EXAMPLE USAGE:
%   buildInternal_mxArrayDef('Internal_mxArray.h')
%
% Author: Bruno Luong <brunoluong@yahoo.com>
%
% History
%   Original: 28-Jun-2009

% Location of the header file matrix.h
MLincludepath = [matlabroot() filesep 'extern' filesep 'include'];
matrixhfile = 'matrix.h';

fid = fopen([MLincludepath filesep matrixhfile]);
if fid>0
    c = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
    try
        fclose(fid);
    end
    
    content = c{1};
    
    % Look for the line containing "struct mxArray_tag {"
    idxmxArray_tag = strfind(content,'struct mxArray_tag {');
    l1 = find(~cellfun('isempty',idxmxArray_tag),1,'first');
    if isempty(l1)
        error('Cannot parse matrix.h file');
    end
    
    % Modify the mxArray_tag to typedef definition
    content{l1} = strrep(content{l1}, ...
                         'struct mxArray_tag', 'typedef struct');
    
    % Loop on the line and stop when the last curly bracket after
    % find the corresponding closed curly bracket
    % "struct mxArray_tag { ... }"
    l9 = 0;
    ncurlybrackets = 0;
    nlevels = 0;
    for l=l1:length(content)
        line  = content{l};
        nopen = sum(line=='{');
        nclose = sum(line=='}');
        ncurlybrackets =  ncurlybrackets + (nopen + nclose);
        nlevels = nlevels + (nopen - nclose);
        if ncurlybrackets>0 && nlevels==0
            l9 = l;
            % Modify the last line with the typedef name 'Internal_mxArray'
            lastcurly = find(line=='}',1,'last');
            line = [line(1:lastcurly) ...
                    ' Internal_mxArray' ...
                    line(lastcurly+1:end)];
            content{l} = line;
            break;
        end
    end
    if l9==0
        error('Cannot parse matrix.h file');
    end
    % Here is the definition we are interested in
    content = content(l1:l9);    

    if nargin>=1
        thisfile = mfilename();
        fid = fopen(mxArraydefFilename,'wt');
        % Write a comment header
        fprintf(fid, ['/* Built automatically by ' thisfile '.m\n']);
        fprintf(fid, ['\tBuilt date: ' datestr(now) '\n']);
        fprintf(fid, ['\tMatlab version: ' version('-release') '\n']);
        fprintf(fid, '*/\n\n');
        
        if fid>0
            % Write the content to header file
            for l=1:length(content)
                fprintf(fid, '%s\n', content{l});
            end
            try
                fclose(fid);
            end
        else
            error('Cannot write the header file %s', mxArraydefFilename);
        end
    end
else % fail to open matrix.h
    error('Cannot find ML <matrix.h> file');
end