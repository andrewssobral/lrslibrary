% function [res loc] = minkmex(list, k)
%
% Matlab C-Mex
% Purpose: Same as MINK, i.e.,
% Return in RES the K smallest elements of LIST
% LOC is Location of the smallest values: RES=LIST(LOC)
% This MEX works on double only, and output RES is unsorted
% Algorithm according to http://en.wikipedia.org/wiki/Selection_algorithm
% Compilation: mex -O -v minkmex.c
% Author Bruno Luong <brunoluong@yahoo.com>
% Last update: 07/April/2009
%

error('Mex file not yet compiled. Action: mex -O -v minkmex.c')