function [X,S]=post_noise_removal(X,D,stdev)

% Post Noise Removal Routine
% Given noise free background X (X is fixed), noise in foreground is removed as follows: 
% S = argmin{|S|_1: |S+X-D|_F<=delta}

[n1,n2]=size(D);
n_min = min(n1,n2);
[U,Sigma_X,V] = svd(X,'econ');
Sigma_X = diag(Sigma_X);
for i=1:n_min-1
    s(i)=Sigma_X(i)/Sigma_X(i+1);
end

[dummy, rr] = max(s);
X=U(:, 1:rr) * diag(Sigma_X(1:rr)) * V(:, 1:rr)';
delta = sqrt(n1*n2)*stdev;
S= l1proj(D-X,delta);