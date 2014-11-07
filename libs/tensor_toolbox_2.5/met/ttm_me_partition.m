function [edims, sdims] = ttm_me_partition(U, esz, n)
%TTM_ME_PARTITION Finds best order for ttm_me.
%
%   [EDIMS, SDIMS] = TTM_ME_PARTITION(U, ESZ, N) finds the best order for
%   memory efficient ttm.  ESZ specifies the number of dimensions that
%   are to be handled elementwise. U is the cell array of matrices.
%   The result returned in EDIMS and SDIMS are the modes that need to
%   be handled in element-wise and in slice-wise respectively. The
%   orders are in the descending order of the reduction ratio.
%
%   See also TUCKER_ME, TTM_ME.
%
%   Based on the paper:
%   T. G. Kolda and J. Sun. Scalable Tensor Decompositions for Multi-aspect
%   Data Mining. In: ICDM 2008: Proceedings of the 8th IEEE International
%   Conference on Data Mining, December 2008.  

% $Id: ttm_me_partition.m,v 1.2 2009/07/08 00:03:52 tgkolda Exp $

if nargin<3
    %default: n is not in 1:length(U)
    n = -1; 
end

% compute reduction ratios
dims = setdiff(1:length(U), n);
r = [];
for i=dims
    [nr, nc] = size(U{i});
    r(end+1) = nr/nc;
end

%sort the reduction ratio and
%put the ones w/ highest reduction ratio in edims
[junk, ind] = sort(r,'descend');
edims = sort(dims(ind(1:min(length(ind),esz))));
sdims = sort(dims(ind(esz+1:end)));
end       