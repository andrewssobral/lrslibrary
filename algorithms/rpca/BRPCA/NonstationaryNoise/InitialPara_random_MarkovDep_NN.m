function  Theta0 = InitialPara_random_MarkovDep_NN(X,K)
if nargin<2
    K=150;
end
[P,N]=size(X);   

%--------Initialize Parameter D,S,Z,Pi,gamma_epsi,gamma_s-----------------
D=zeros(P,K);
for k=1:K
    D(:,k)=randn(P,1)*sqrt(1/P);
end


S=zeros(K,N);
for n=1:N
    S(:,n)=randn(K,1);
end


Z = zeros(K,1);
%Z = ones(K,1);
Delta = ones(K,1);
Tao = ones(K,1);
Pi = 0.5*ones(K,1);         
gamma_epsi = 1e3*ones(1,N);  

%sampe sparse component
gamma_s = 1;
S2 = randn(P,N);
Z2 = zeros(P,N);
Pi2 = 0.5*ones(P,N);




Theta0.D = D;
Theta0.S = S;
Theta0.Z = Z;
Theta0.Delta = Delta;
Theta0.Tao = Tao;
Theta0.Pi = Pi;
Theta0.gamma_epsi = gamma_epsi;
Theta0.S2 = S2;
Theta0.Z2 = Z2;
Theta0.gamma_s = gamma_s;
Theta0.Pi2 = Pi2;
end