%this is the code for AAAI 2016 paper: Top-N Recommender System via Matrix
%Completion by zhao kang,chong peng,qiang cheng.

%[hr,arhr,X] = mc_logdet(data,mu,rho)
function [X] = demo(data,mu,rho,toler,maxiter)
  [m,n]=size(data);
  M=data;
  id=find(M==0);
  ID=ones(m,n);
  ID(id)=0;
  [X] = MC_LogDet_v3(M,ID,mu,rho,toler,maxiter);
end




