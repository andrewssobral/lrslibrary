% Serhat Selcuk Bucak, bucakser@msu.edu

% Lets assume that all samples are stored in data matrix V (each column is a sample)

%Execute NMF for the first n samples
%% These number are selected just for demonstration
V=rand(40,500);
n=100;
rdim=10;
maxiter=150;
[W, H, objhistory] = nmf(V(:,1:n), rdim, 0, maxiter);
% Now we can execute inmf on each new samples
maxiter=50;
A=V(:,1:n)*H';
B=H*H';
h=H(:,end); % Warm start for h
for i=n+1:size(V,2)
    i
    V_new=V(:,i);
    [W_new, h, A, B] = inmf( V_new, W, h, A, B, rdim, 0.9, 0.1, maxiter);
    H_store(:,i-n)=h; %Just for demonstration
end