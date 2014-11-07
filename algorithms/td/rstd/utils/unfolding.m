% unfolding the tensor along the kth dimension
% use permute for efficiency, this is much faster than circshift
function [matrix] = unfolding(tensor, k)

site = size(tensor);

if k==1
    matrix = reshape(tensor, site(1), []);
    return;
end

dim = 1:ndims(tensor);
dim = dim';

dim = circshift(dim, 1-k);
matrix = reshape(permute(tensor, dim), site(k), []);



