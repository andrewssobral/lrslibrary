function U = DLA(fea, gnd, k1, k2, beta)

% Discriminative Locality Alignment (DLA)
% Reference:
%   T. Zhang, D. Tao, X. Li, and J. Yang, "Patch Alignment for
%   Dimensionality Reduction," IEEE Transactions on Knowledge and Data Engineering, vol. 21, no. 9, pp. 1299–1313, 2009.
%
% Arguments:
% <Input>:
%   fea (nSmp x nFea): observation matrix;
%   gnd: label for each observation;
%   k1: nearest neighbor number on within-class patch;
%   k2: nearest neighbor number on between-class patch;
%   beta: trade-off between two patches.
% <Output>:
%   U (nFea x nFea): whole mapping matrix, where U(:,1:r) selects first r
%   atoms to construct the mapping function.
%
% Written by Naiyang Guan (ny.guan@gmail.com)

if k2 <= 0
    fprintf('error! k2 must be positive.\n');
    return;
end
N = size(fea, 1);
dist = EuDist2(fea,fea,0);

L = zeros(N);
omega = [ones(k1,1); -beta*ones(k2,1)];
Li = [sum(omega), -omega'; -omega, diag(omega)];

for i = 1 : N
    sam_idx = find(gnd == gnd(i));
    diff_idx = find(gnd ~= gnd(i));
    sam_dist = dist(i, sam_idx');
    diff_dist = dist(i, diff_idx');
    [junk, index] = sort(sam_dist);
    sam_idx = sam_idx(index);
    [junk, index] = sort(diff_dist);
    diff_idx = diff_idx(index);
    
    Fi = [sam_idx(1:(k1+1)); diff_idx(1:k2)];
    
    L(Fi, Fi) = L(Fi, Fi) + Li;
end

data = fea'*L*fea;
data = max(data, data');

[U, eig_value] = eig(data);
eig_value = diag(eig_value);
[junk, index] = sort(eig_value);
U = U(:, index);

% normalization
U = U*diag(1./((sum(U.^2).^0.5)));

return;