function Y=rtr_l0(X,s)

[n1,n2]=size(X);
X=X(:);
i1=find(X~=0);
Y=X(i1);
[~,i2]=sort(abs(Y(:)),'descend');
X(i1(i2(s+1:end)))=0;
Y=reshape(X,n1,n2);


%%%% old version %%%
% [n1,n2]=size(X);
% 
% Y=X(:);
% [~,ii]=sort(abs(Y(:)));
% Y(ii(1:end-s))=0;
% Y=reshape(Y,n1,n2);

