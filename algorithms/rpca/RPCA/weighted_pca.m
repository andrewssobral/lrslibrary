function [Bg,cg,info,meanw]=missing_bc(Data,W,components,max_iter,display,Bg,cg,meanw)
%function [Bg,cg,info,meanw]=missing_bc(Data,W,components,max_iter,display,Bg,cg,meanw)
%Usage: 
%	Input:
%			Data: Data matrix.
%			Components: Number of components to extract.
%			max_iter: Maximum number of iterations.
%			diplay_iter: Number of iterations to display the information.
%			Bg: Initial estimation of the bases.
%			cg: Initial estimation of the coefficients.
%			meanw: Initial estimation of the mean.
%	Output:
%			Bg: Final estimation of the bases.
%			cg: Final estimation of the coefficients.
%			info: Information about the error.
%			meanw: Final estimation of the mean.
% Copyright 2001 Fernando De la Torre 2001 


%inicialize variables
info=zeros(max_iter,2);

iter=1;flag=1;

while  iter<max_iter & flag
   
error=(Data-meanw*ones(1,size(Data,2))-Bg*cg);
errsvd=sum(sum((error.^2).*W));

clear error

cgant=cg; Bgant=Bg;

%computing the weighted mean
meanw= sum((W.*(Data-Bg*cg))')'./sum(W')';

%computing the weighted coefficients
for i=1:size(Data,2)
      Ap=W(:,i)*ones(1,size(Bg,2)).*Bg;
      cg(:,i)=inv(Bg'*Ap)*Ap'*(Data(:,i)-meanw);
end;

%computing the weighted bases
for i=1:size(Bg,1)
   Ap=cg.*(ones(size(cg,1),1)*W(i,:));
   Bg(i,:)=(inv(cg*Ap')*Ap*(Data(i,:)-meanw(i)*ones(1,size(Data,2)))')';   
end

%computing the angular error
angular_error=subspace(Bg,Bgant);
if angular_error<1e-3 &(iter>30)
      	flag=0;
end
    
if  rem(iter,display)==0
   fprintf('Iter:%d , Err:%.3f ,  angular_error: %.3f \n',iter,errsvd, angular_error);   
end;

info(iter,:)=[errsvd angular_error];
iter=iter+1;

end
return;

