%% initialization

[U,S,V]=svd(Z,'econ');
% [U,S,V]=svd(randn(n1,n2),'econ');
U=U(:,1:n3);
V=V(:,1:n3);
S=S(1:n3,1:n3);
        
A=U*S*V';
B=rtr_l0(Z-A,n4);
