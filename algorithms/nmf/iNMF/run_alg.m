% iNMF: Incremental Subspace Learning via NMF (Bucak and Gunsel, 2009)
% process_video('NMF', 'iNMF', 'dataset/demo.avi', 'output/demo_iNMF.avi');

% Execute NMF for the first n samples
n = 10; % first 10 samples
rdim = 1; % rank-1
maxiter = 150;
[W,H] = nmf(M(:,1:n), rdim, 0, maxiter);
L = W*H;
% Now we can execute iNMF on each new samples
maxiter = 50;
A = M(:,1:n)*H';
B = H*H';
h = H(:,end); % Warm start for h
for i = n+1:size(M,2)
  disp(i);
  M_new = M(:,i);
  [W_new,h,A,B] = inmf(M_new,W,h,A,B,rdim,0.9,0.1,maxiter);
  % H_store(:,i-n) = h; % Just for demonstration
  L(:,end+1) = W_new*h;
end
S = M - L;