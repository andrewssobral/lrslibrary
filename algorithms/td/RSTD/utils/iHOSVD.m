function [data] = iHOSVD(core, U)

site = size(core);
data = core;

for i=1:ndims(core)
    [m n] = size(U{i});
    site(i) = m;
    data = folding(U{i}*unfolding(data,i), i, site);
end