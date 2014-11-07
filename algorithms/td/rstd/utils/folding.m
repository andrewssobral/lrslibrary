% folding the matrix into a tensor along the kth dimension
% use permute for efficiency, this is much faster than circshift

function [tensor] = folding(matrix, k, site)

dim = 1:length(site);
dim = dim';

dim = circshift(dim, k-1);
site = (circshift(site', 1-k))';

tensor = permute(reshape(matrix, site), dim);