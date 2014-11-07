function [core] = HOSVD(data, U)

site = size(data);
core = data;

for i=1:ndims(data)
    [m n] = size(U{i});
    site(i) = n;
    core = folding(U{i}'*unfolding(core,i), i, site);
end
