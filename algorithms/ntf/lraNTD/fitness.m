function [fit res]=fitness(Y,Ycap,flag)
%% function [fit res]=fitness(Y,Ycap)

if strcmp(class(Y),'double')
    Y=tensor(Y);
end


normY=norm(Y);
normYcap=norm(Ycap);
res = abs(sqrt( normY^2 + normYcap^2 - 2 * innerprod(Y,Ycap) ));
fit=1-res/normY;

if nargin==3 && flag==2
    fit(2)=1-(res/normY).^2;
end