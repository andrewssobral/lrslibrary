function varargout = mink(varargin)
% function res = MINK(list, k)
% 
% If LIST is a vector, MINK returns in RES the K smallest elements of LIST
%   RES is sorted in ascending order
% [res loc] = MINK(...)
%   Location of the smallest elements: RES=LIST(LOC)
% If list is a matrix, MINK operates along the first dimension
% Use MINK(..., dim) to operate along the dimension dim
%     MINK(..., dim, 'sorting', false) to disable the post-sorting step
%                                      (true by default)
%
% Author Bruno Luong <brunoluong@yahoo.com>
% Contributor: Matt Fig
% Last update: 07/April/2009
%              10/Jan/2010: possibility to disable post-sorting step


nout=cell(1,max(1,nargout));
[nout{:}] = minmaxk(@minkmex, varargin{:});
varargout=nout;
