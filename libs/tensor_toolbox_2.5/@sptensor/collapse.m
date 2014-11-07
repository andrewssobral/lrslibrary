function s = collapse(t,dims,fun)
%COLLAPSE Collapse sparse tensor along specified dimensions.
%
%  S = COLLAPSE(T,DIMS) sums the entries of T along all dimensions
%  specified in DIMS. If DIMS is negative, then T is summed across
%  all dimensions *not* specified by -DIMS.
%
%  S = COLLAPSE(T) is shorthand for S = COLLAPSE(T,1:ndims(T)).
%
%  S = COLLAPSE(T,DIMS,FUN) accumulates the entries of T using the
%  accumulation function @FUN.
%
%  Examples
%  subs = [1 1 1; 1 1 3; 2 2 4; 4 4 4]
%  vals = [10.5; 1.5; 2.5; 3.5]
%  X = sptensor(subs,vals,[4 4 4]);
%  Y = collapse(X,[2 3]) %<-- sum of entries in each mode-1 slice
%  Y = collapse(ones(X),[1 2]) %<-- nnz in each mode-3 slide
%  Y = collapse(ones(X),[1 2],@max) %<-- 1 if mode-3 has any entry
%  Y = collapse(ones(X),-3,@max); %<-- equivalent
%
%  See also SPTENSOR, SPTENSOR/SCALE.
%
%MATLAB Tensor Toolbox.
%Copyright 2012, Sandia Corporation.

% This is the MATLAB Tensor Toolbox by T. Kolda, B. Bader, and others.
% http://www.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2012) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in the file LICENSE.txt


if ~exist('fun', 'var')
    fun = @sum;
end

if ~exist('dims', 'var')
    dims = 1:ndims(t);
end

dims = tt_dimscheck(dims,ndims(t));
remdims = setdiff(1:ndims(t),dims);

% Check for the case where we accumulate over *all* dimensions
if isempty(remdims)
    s = fun(t.vals);
    return;
end

% Calculate the size of the result
newsiz = size(t,remdims);

% Check for the case where the result is just a dense vector
if numel(remdims) == 1
    if ~isempty(t.subs)
        s = accumarray(t.subs(:,remdims), t.vals, [newsiz 1], fun);
    else
        s = zeros(newsiz,1);
    end
    return;
end

% Create the result
if ~isempty(t.subs)
  s = sptensor(t.subs(:,remdims), t.vals, newsiz, fun);
else
  s = sptensor([],[],newsiz);
end
