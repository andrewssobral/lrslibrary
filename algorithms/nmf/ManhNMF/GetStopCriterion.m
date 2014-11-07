function retVal=GetStopCriterion(stop_rule,X,gradX)
% Stopping Criterions
% Written by Naiyang (ny.guan@gmail.com)

switch stop_rule
    case 1
        pGrad=gradX(gradX<0|X>0);
        retVal=norm(pGrad);
    case 2
        pGrad=gradX(gradX<0|X>0);
        pGradNorm=norm(pGrad);
        retVal=pGradNorm/length(pGrad);
    case 3
        resmat=min(X,gradX); resvec=resmat(:);
        deltao=norm(resvec,1);  %L1-norm
        num_notconv=length(find(abs(resvec)>0));
        retVal=deltao/num_notconv;
end