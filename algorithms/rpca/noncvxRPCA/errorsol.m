function [error]= errorsol(Y1,X,W,lambda,mu,type)
switch type
  case 1
    D = -Y1/mu+(X-W);
    E = zeros(size(D));
    epsilon = lambda/mu;
    DD = abs(D)-epsilon;
    DD2 = DD.*sign(D);
    ID = abs(D)>epsilon;
    E(ID) = DD2(ID);
  case 21
    alpha=lambda/mu;
    G=X-W-Y1/mu;
    G1 = sqrt(sum(G.^2,1));
    G1(G1==0) = alpha;
    G2 = (G1-alpha)./G1;
    E = G*diag((G1>alpha).*G2);
end
error=E;