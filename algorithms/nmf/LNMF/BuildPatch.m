function [D,W,Lw,Lb]=BuildPatch(fea,gnd,k1,k2,zeta)

% Build two patches for MD-NMF and NDLA by using KNN + sparse binary model
%   data: nSample x nFeature matrix,
%   zeta: perturbution coefficient.
% Written by Naiyang Guan (ny.guan@gmail.com)
% Copyright by Naiyang Guan and Dacheng Tao

% Initialization: assume data is MxN matrix
N=size(fea, 1);
Sw=zeros(N);
Sb=zeros(N);

% Extreme position
if k2<=0
    fprintf('error! k2 cannot be smaller than 1\n');
    return;
end

% Euclidean distance between sample points in the PCA space
dist=EuDist2(fea,fea,0);

% Construct sparse Laplacian matrix on (k1,k2)-patch
for i=1:N
    % k1 nearest neighbors within same class
    sam_idx=find(gnd==gnd(i));
    sam_dist=dist(i,sam_idx');
    [junk,index]=sort(sam_dist);
    sam_idx=sam_idx(index);
    if k1>=1
        Fi=sam_idx(2:(k1+1));
    else
        Fi=i;
    end
    Sw(i,Fi)=1;
    Sw(Fi,i)=1;
    
    % k2 nearest neighbors between classes
    diff_idx=find(gnd~=gnd(i));
    diff_dist=dist(i,diff_idx');
    [junk,index]=sort(diff_dist);
    diff_idx=diff_idx(index);
    Fi=diff_idx(1:k2);
    Sb(i,Fi)=1;
    Sb(Fi,i)=1;
end

% Laplacian regularization and normalization
Dw=diag(sum(Sw));
Db=diag(sum(Sb));
Lw=Dw-Sw;
Lb=Db-Sb+eye(N)*sum(sum(Sb))*zeta;

[U,S,V]=svd(Lb);
S=diag(S);
S=diag(1./sqrt(S));
D=U*S*V'*Dw*V*S*U';
W=U*S*V'*Sw*V*S*U';

return;