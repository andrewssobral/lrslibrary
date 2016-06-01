function [I J col omega] = my_randsample(m, n, k)
%generates uniformly distribured samples at k positions of an mxn matrix

omega = ceil(rand(k, 1) * m * n);
omega = unique(omega);
while length(omega) < k    
    omega = [omega; ceil(rand(k-length(omega), 1)*m*n);];
    omega = unique(omega);
end

[I,J] = ind2sub([m,n], omega);
col = [0; find(diff(J)); k]; 




