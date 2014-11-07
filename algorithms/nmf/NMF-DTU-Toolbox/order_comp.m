function [W,H,nrgy]=order_comp(W,H)
%
% Order components according to "energy"
%
[D,K]=size(W);
[K,N]=size(H);
%
nrgy=zeros(K,1);
wsum=sum(W,1);
hsum=sum(H,2)';
nrgy=wsum.*hsum;
[nrgy,index]=sort(-nrgy);
nrgy=-nrgy;
W=W(:,index);
H=H(index,:);