function [max_mem] = ttm_me_mem(X, U, edims, sdims, tflag)
%TTM_ME_MEM Estimates intermediate memory comsumption for ttm_me.
%
%  MAX_MEM = TTM_ME_MEM(X, U, EDIMS, SDIMS, TFLAG) estimates the memory
%  comsumption used in the TTM_ME function. Here, X is a sparse tensor  
%  (sptensor), U is a cell array of matrices of length ndims(X),
%  EDIMS specifies the dimensions that are to be handled elementwise, sdims
%  specifies those dimensions that are to be handled in the standard way,
%  and TFLAG indicates to multiply the matrix transpose.
%
%   See also TENSOR_TOOLBOX, TUCKER_ME, TTM_ME.
%
%   Code by Tamara Kolda and Jimeng Sun, 2008. 
%
%   Based on the paper:
%   T. G. Kolda and J. Sun. Scalable Tensor Decompositions for Multi-aspect
%   Data Mining. In: ICDM 2008: Proceedings of the 8th IEEE International
%   Conference on Data Mining, December 2008.  

% $Id: ttm_me_mem.m,v 1.2 2009/07/08 00:03:52 tgkolda Exp $

%% Setup and error checking

% Check number of input arguments
if (nargin < 4)
    error('TTM_ME requires four arguments');
end
if (nargin == 4)
    tflag = '';
else
    tflag = 't';
end
% Check that X is a sparse tensor
if ~isa(X, 'sptensor')
    error('Input tensor X must be sparse');
end

% Set number of dimensions of X
N = ndims(X);

% Check that U is a cell array
if ~iscell(U)
    error('U must be a cell array');
end

% Check that the cell array U
if numel(U) ~= N
    error('Incorrect number of elements in U');
end

% Check that every member of edims is in 1:N
tf = ismember(edims, 1:N);
if min(tf) == 0
    error('Invalid dimensions specified');
end

% Check that every member of sdims is in 1:N
tf = ismember(sdims, 1:N);
if min(tf) == 0
    error('Invalid dimensions specified');
end

% Check that edims and sdims are distinct
idx = intersect(edims, sdims);
if ~isempty(idx)
    error('Invalid dimensions specified');
end

%% Calculate the max memory from two sources
% 1. intermediate tensor
    sizY = size(X);    
    for n = edims
        sizY(n) = 1;    
    end
    % special case for the MET(0) or original Tucker-ALS
    if isempty(edims)
        if strcmp(tflag, 't')
            sizY(sdims(1)) = size(U{sdims(1)}, 2);
        else
            sizY(sdims(1)) = size(U{sdims(1)}, 1);
        end
    end
    mem = prod(sizY);    
% 2. ttv 
    if ~isempty(edims)
        mem = mem + nnz(X);
    end
    max_mem = mem;
end