function grasta_path(BaseDirectory)
    fprintf('Setting necessary paths ... \n');

    if nargin==0
        Prefix  = [pwd filesep];
    else
        Prefix  = [BaseDirectory filesep];
    end    
    appendpath(Prefix);
    appendpath([Prefix 'Mex']);
    appendpath([Prefix 'util']);
    appendpath([Prefix 'make_video']);
    
    fprintf('Disabling case sensitivity warning ... \n');
    warning('off','MATLAB:dispatcher:InexactMatch');
end

function appendpath(string)
    fprintf('\t%s \n', string);
    addpath(genpath(string));
end


