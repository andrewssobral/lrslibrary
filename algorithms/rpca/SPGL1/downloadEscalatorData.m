function downloadEscalatorData
if ~exist('escalator_data.mat','file')
    % This downloads a 3.7 MB demo
    disp('About to download the 3.7 MB video...');
    
    baseDirectory = fileparts(mfilename('fullpath'));
    disp(baseDirectory)
    
    urlwrite('http://cvxr.com/tfocs/demos/rpca/escalator_data.mat',...
        fullfile(baseDirectory,'escalator_data.mat'));
else
    disp('You already have the data file; everything is good');
end