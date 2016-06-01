% generate path
path(path,genpath(pwd))

% compile mex files if necessary
Amexfile = ['partXY.mex' lower(computer)];
if ~exist(Amexfile,'file')
    cd Utilities;
    %mex -O fastWHtrans.cpp
    %mex -O -largeArrayDims updateSval.c -lmwblas -output  updateSval
    mex -O -largeArrayDims updateSval.c 
    %mex -O updateSvalZw.c
    %mex -O -largeArrayDims partXY.c -lmwblas -output partXY
    mex -O -largeArrayDims partXY.c
    %mex -O proj_LR_sum.c
    cd ..;
end

% you are all set
disp('Welcome to LMaFit')
