function [X,S]=treshold(X,D,stdev)
%Treshold routine
[n1,n2]=size(D);
n_min = min(n1,n2);
[U,Sigma_X,V] = svd(X,'econ');
Sigma_X = diag(Sigma_X);
for i=1:n_min-1
    s(i)=Sigma_X(i)/Sigma_X(i+1);
end

[dummy, rr] = max(s);
X = U(:, 1:rr) * diag(Sigma_X(1:rr)) * V(:, 1:rr)';
%Sigma_X(rr+1:n_min)=0;
deltabar = sqrt(n1*n2)*stdev;
S = l1proj(D-X,deltabar);

% rel_err_X = norm(X-optimal_X,'fro')/norm(optimal_X,'fro')
% Sigma_opt = svd(optimal_X, 'econ');
% ind_pos_X = find(Sigma_opt>1e-10);
% ind_0_X = find(Sigma_opt<=1e-10);
% inf_err_pos_X = norm(Sigma_X(ind_pos_X)-Sigma_opt(ind_pos_X),inf);
% inf_err_0_X = norm(Sigma_X(ind_0_X)-Sigma_opt(ind_0_X),inf);
% 
% display(['max{|sigmaX(i)-sigmaX*(i)|: i s.t. sigmaX*(i)>0}: ',num2str(inf_err_pos_X)])
% display(['max{|sigmaX(i)-sigmaX*(i)|: i s.t. sigmaX*(i)=0}: ',num2str(inf_err_0_X)])
% 
% rel_err_S = norm(S-optimal_S,'fro')/norm(optimal_S,'fro')
% ind_pos_S = find(abs(optimal_S)>1e-10);
% ind_0_S = find(abs(optimal_S)<=1e-10);
% inf_err_pos_S = norm(S(ind_pos_S)-optimal_S(ind_pos_S),inf);
% inf_err_0_S = norm(S(ind_0_S)-optimal_S(ind_0_S),inf);
% 
% display(['max{|S(ij)-S*(ij)|: ij s.t. S*(ij)>0}: ',num2str(inf_err_pos_S)])
% display(['max{|S(ij)-S*(ij)|: ij s.t. S*(ij)=0}: ',num2str(inf_err_0_S)])