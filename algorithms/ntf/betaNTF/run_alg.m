% NTF | betaNTF | Simple beta-NTF implementation (Antoine Liutkus, 2012)
% process_video('NTF', 'betaNTF', 'dataset/demo.avi', 'output/demo_beta-NTF.avi');

% Compute a simple NTF model of 10 components
A = double(T);
r = 10;
[~,~,~,L] = betaNTF(A,r);
S = (A - L);
% For reconstruction
% for i = 1:r, B_hat(:,:,i) = W * diag(Q(i,:)) * H'; end
