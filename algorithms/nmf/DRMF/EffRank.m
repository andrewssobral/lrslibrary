function [ r ] = EffRank(s, thresh)
%[ r ] = EffRank(s, thresh)
% get the effective rank of an singular value array
% s: the singular values in descending order
% thresh: the total variance to preserve
% author: Liang Xiong (lxiong@cs.cmu.edu)

if nargin < 2; thresh = 0.99; end
assert(all(s > 0), 'not a proper singular value array');

s = sort(s.^2,'descend');
cs = cumsum(s);
if cs(end) < 1e-10
    error('EffRank has encounted a zero matrix');
else
    r = sum(cs./cs(end) <= thresh) + 1;
    r = min(length(s), r);
end
