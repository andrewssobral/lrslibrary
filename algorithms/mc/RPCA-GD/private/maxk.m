function varargout = maxk(varargin)
% function res = MAXK(list, k)
% 
% If LIST is a vector, MAXK returns in RES the K largest elements of LIST
%   RES is sorted in descending order
% [res loc] = MAXK(...)
%   Location of the largest elements: RES=LIST(LOC)
% If list is a matrix, MAXK operates along the first dimension
% Use MAXK(..., dim) to operate along the dimension dim
%     MAXK(..., dim, 'sorting', false) to disable the post-sorting step
%                                      (true by default)
%
% Author Bruno Luong <brunoluong@yahoo.com>
% Contributor: Matt Fig
% Last update: 07/April/2009
%              10/Jan/2010: possibility to disable post-sorting step

nout=cell(1,max(1,nargout));
[nout{:}] = minmaxk(@maxkmex, varargin{:});
varargout=nout;

