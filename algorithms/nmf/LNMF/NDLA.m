function [W,H,ela]=NDLA(fea,gnd,k,r,beta,verbose)

% Non-negative Discriminative Locality Alignment (NDLA)
% Written by Naiyang Guan (ny.guan@gmail.com)
% Copyright 2012 by Naiyang Guan and Dacheng Tao

% Reference:
% [1] N. Guan, D. Tao, Z. Luo, and B. Yuan, "Manifold Regularized Discriminative Non-negative Matrix Factorization with Fast Gradient Descent,"
% IEEE Transactions on Image Processing, vol. 20, no. 7, pp. 2030–2048, 2011.
% [2] N. Guan, D. Tao, Z. Luo, and B. Yuan, "Non-negative Patch Alignment Framework,"
% IEEE Transactions on Neural Networks, vol. 22, no. 8, pp. 1218–1230, 2011.

% Arguments:
% <Input>:
%   fea: data matrix (nSample x nFeature);
%   gnd: label information;
%   k: number of nearest neighbors in within-class and between-class
%   patches;
%   r: lower dimensionality;
%   beta: trade-off parameter,
% <Output>:
%   W: basis matrix;
%   H: coefficient matrix.
% <Usage Example>:
%   >>[W,H]=NDLA(V,gnd,[2,15],50,1e-3);

% Build Patches on Within-class and Between-class Graphs
fea=fea/max(fea(:));
zeta=1e-4;      % Factor for positive definite perturbution
[D,S]=BuildPatch(fea,gnd,k(1),k(2),zeta);
[nSmp,nFea]=size(fea);
V=fea';

% Initialization by One-step Normal MUR
W0=rand(nFea,r);
H0=rand(r,nSmp);
H0=H0.*(H0*S+W0'*(V./(W0*H0)))./(H0*D+sum(W0)'*ones(1,nSmp));
W0=W0.*((V./(W0*H0))*H0')./(ones(nFea,1)*sum(H0,2)');

% NDLA by Required Algorithm
[W,H,it,ela,HIS]=NPAF(V,D,S,r,'beta',beta,'alg_type','smur','w_init',W0,'h_init',H0,'verbose',verbose);

% Report!
if verbose,
    fprintf('NDLA success with the following details:\n');
    fprintf('\titeration number: %d,\n',it);
    fprintf('\telapse time: %.3f,\n',ela);
    fprintf('\tfinal objective value: %.3e.\n',HIS.objs(end));
end

return;