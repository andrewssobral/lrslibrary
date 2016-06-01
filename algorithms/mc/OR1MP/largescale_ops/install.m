clear; clc;
disp('compile native libraries')
mex -O -largeArrayDims sparse_inp.c
mex -O -largeArrayDims sparse_update.c
disp('done compiling.')
test_consist1();
test_consist2();