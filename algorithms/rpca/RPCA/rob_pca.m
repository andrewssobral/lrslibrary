function [rob_mean,Bg,cg,info,Sigma]=rob_bc_nw(Data,components,max_iter,display,Sigma_end,Sigma_start,iter_grad,Bg,cg,rob_mean)
%function [rob_mean,Bg,cg,info,Sigma]=rob_bc_nw(Data,components,max_iter,display,Sigma_end,Sigma_start,iter_grad,Bg,cg,rob_mean)
%Usage: 
%	Input:
%			Data: Data matrix (each column is a vectorized image).
%			Components: Number of components to extract.
%			max_iter: Maximum number of iterations.
%			diplay_iter: Display results every display_iter iterations.
%			Sigma_end: Final sigma.
%			Sigma_start: Initial sigma.
%			iter_grad: Iterations of the internal update equations. 
%	Output:
%			rob_mean: Robust mean.
%			Bg: Columns: Robust bases.
%			cg: Columns: Robust coefficients.
%			info: Matrix containing : [ robust_error angular_error time flops].
%			 
% Robust Principal Component Analysis/ Singular Value Decomposition 
% Copyright 2001 Fernando De la Torre 2001.
 


%inicialize variables
%info=zeros(max_iter,4);
info=zeros(max_iter,2);
if nargin==7
   rob_mean=mean(Data')';
   [u,s,v]=svd(Data-rob_mean*ones(1,size(Data,2)),0);
   Bg=u(:,1:components);
   cg=(Bg'*(Data-rob_mean*ones(1,size(Data,2))));
   clear u s v  
end

Sigma=Sigma_start; %initial Sigma value 
	
iter=1; 
flag=1;mu=1;
errsvd=0;

while  iter<max_iter & flag
   
sigma_max_ant=max(Sigma);
sigma_min_ant=min(Sigma);

if (rem(iter,1)==0)
	Sigma=Sigma*0.92;   
   Sigma=max(Sigma,Sigma_end);
end
   
Sigmat=Sigma.^2*ones(1,size(Data,2));

%flops(0)
%tic

cgant=cg; Bgant=Bg;
%updating the mean
for i=1:iter_grad
   error=(Data-rob_mean*ones(1,size(Data,2))-Bg*cg);
   pondera=(error.*Sigmat)./((Sigmat+error.^2).^2);   
   rob_mean=rob_mean+mu*sum(pondera')'./(size(Data,2)*1./Sigmat(:,1));
end  

%computing the error
errsvdant=errsvd;
errsvd=sum(sum(((error.^2)./(Sigmat+error.^2)) ));  

%updating the basis
for i=1:iter_grad
   error=(Data-rob_mean*ones(1,size(Data,2))-Bg*cg);
   pondera=error.*((Sigmat))./((Sigmat+error.^2).^2);   
   Bg=Bg+mu*(pondera*cg')./((1./Sigmat)*(cg.*cg)');
end

%updating the coefficients
for i=1:iter_grad
  error=(Data-rob_mean*ones(1,size(Data,2))-Bg*cg);
  pondera=error.*((Sigmat))./((Sigmat+error.^2).^2);   
  cg=cg+mu*(Bg'*pondera)./((Bg.*Bg)'*(1./Sigmat));
end

%time=toc;
%compute=flops;
   
 
angular_error=subspace(Bg,Bgant);
    
if  rem(iter,display)==0
   fprintf('Iter:%d , Err:%.3f , Sigmax:%.4e, Sigmin:%.4e,  angular_error: %.4f mu %.2f \n',iter,errsvd,max(Sigma),min(Sigma), angular_error,mu);
end;   
    
%info(iter,:)=[errsvd angular_error time compute];
info(iter,:)=[errsvd angular_error];

if (errsvd>errsvdant) & (sigma_max_ant==max(Sigma)) & (sigma_min_ant==min(Sigma))
   mu=mu*0.9;
end

if angular_error<1e-4 & (iter>30)
      	flag=0;
end
      
iter=iter+1;
end;

return;

