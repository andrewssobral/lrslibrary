function [in_perm,error]=compare(W1,W2)
%
% This functions compares two W'matrices from NMF
% by estimating the permutation and computing the
% normalized LS of the permuted matrix
%
% in_perm is the permutation of W2 to similar it with W1
% error is 1-mean(correff);
[D,K]=size(W1);
covar=zeros(K);
xW1=W1-repmat(mean(W1,1),D,1);
xW2=W2-repmat(mean(W2,1),D,1);
%
stds1=std(W1,[],1);
stds2=std(W2,[],1);
%    
for k1=1:K,
    %
    for k2=1:K,
        covar(k1,k2)=abs(sum((xW1(:,k1).*xW2(:,k2)))/D)/(stds1(k1)*stds2(k2));
   end,
end,
%
error=0;
inlist=logical(ones(K,1));
arrayx=1:K;
for k=1:K,
   [dummy,in]=max(covar(k,inlist));
   xarrayx=arrayx(inlist);
   in_perm(k)=xarrayx(in);
   inlist(xarrayx(in))=0;
   error=error+dummy;
end
error=1-(error/K);
   