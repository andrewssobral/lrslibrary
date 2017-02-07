% function [res loc] = maxkmex(list, k)
%
% Matlab C-Mex
% Purpose: Same as MAXK, i.e.,
% Return in RES the K largest elements of LIST
% LOC is Location of the largest values: RES=LIST(LOC)
% This MEX works on double only, and output RES is unsorted
% Algorithm according to http://en.wikipedia.org/wiki/Selection_algorithm
% Compilation: mex -O -v maxkmex.c
% Author Bruno Luong <brunoluong@yahoo.com>
% Last update: 07/April/2009
%

error('Mex file not yet compiled. Action: mex -O -v maxkmex.c');