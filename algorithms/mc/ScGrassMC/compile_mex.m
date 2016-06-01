fprintf('Compiling mex functions\n');
mex -largeArrayDims maskmult.c
mex -largeArrayDims setsparseentries.c
