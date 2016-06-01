% generate path
path(path,genpath(pwd))

% compile mex files if necessary
Amexfile = ['partXY.mex' lower(computer)];
if ~exist(Amexfile,'file')
    mex -O -largeArrayDims updateSval.c 
    mex -O -largeArrayDims partXY.c
end

% you are all set
disp('Welcome to mcnf')
