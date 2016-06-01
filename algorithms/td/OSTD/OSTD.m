%{
Online Stochastic Tensor Decomposition for Background Subtraction

Author: Andrews Sobral
https://github.com/andrewssobral/ostd

Reference:
@inproceedings{ostd,
author    = {Sobral, Andrews and Javed, Sajid and Ki Jung, Soon and Bouwmans, Thierry and Zahzah, El-hadi},
title     = {Online Tensor Decomposition for Background Subtraction in Multispectral Video Sequences},
booktitle = {IEEE International Conference on Computer Vision Workshops (ICCVW)},
address   = {Santiago, Chile},
year      = {2015},
month     = {December},
url       = {https://github.com/andrewssobral/ostd}
}
%}
function [Tlowrank,Tsparse,Tmask,Tm] = OSTD(T,k,Tm)
  Tlowrank = zeros(size(T));
  Tsparse = zeros(size(T));
  nmode = length(size(T));
  for j = 1:nmode
    %disp(['tensor.mode-' num2str(j)]);
    Titm = tenmat(T,j); M = double(Titm); % imagesc(Ti)
    z = M(:); %disp(size(z))
    lambda1 = (1/sqrt(max(size(z))));
    lambda2 = lambda1*10;
    if(k == 1)
      [ndim,~] = size(z);
      nrank = 10;
      Tm{j}.A = zeros(nrank,nrank);
      Tm{j}.B = zeros(ndim,nrank);
      
      %%% Initialization by Random Basis
      %Tm{j}.L = rand(ndim,nrank);
      
      %%% Initialization by Bilateral Random Projections
      power = 1;
      Y2 = randn(size(z,2),nrank);
      for jj = 1:power+1
        Y1 = z*Y2;
        Y2 = z'*Y1;
      end
      [Q,R] = qr(Y2,0);
      L = (z*Q)*Q';
      L = repmat(L,[1 nrank]);
      Tm{j}.L = L;
    end
    A = Tm{j}.A;
    B = Tm{j}.B;
    L = Tm{j}.L;
    [r,e] = solve_proj2(z, L, lambda1, lambda2);
    A = A + r * r';
    B = B + (z-e) * r';
    L = update_col_orpca(L, A, B, lambda1);
    lowrank = reshape(L*r,size(M,1),size(M,2));
    sparse = reshape(e,size(M,1),size(M,2));
    Ti_lowrank = double(tensor(tenmat(lowrank,Titm.rdims,Titm.cdims,Titm.tsize)));
    Ti_sparse = double(tensor(tenmat(sparse,Titm.rdims,Titm.cdims,Titm.tsize)));
    Tm{j}.A = A;
    Tm{j}.B = B;
    Tm{j}.L = L;
    Tlowrank = Tlowrank + Ti_lowrank;
    Tsparse = Tsparse + Ti_sparse;
  end % end tensor mode-i
  Tlowrank = Tlowrank./nmode;
  Tsparse = Tsparse./nmode;
  Tmask = medfilt2(double(hard_threshold(mean(Tsparse,3))),[5 5]);
end
