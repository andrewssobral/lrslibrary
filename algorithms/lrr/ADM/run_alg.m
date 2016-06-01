% LRR | ADM | Alternating Direction Method (Lin et al. 2011)
% process_video('LRR', 'ADM', 'dataset/demo.avi', 'output/demo_LRR-ADM.avi');

lambda = 0.1;
rho = 1.9;
DEBUG = 1;
% X = XZ+E
[Z,E] = adm_lrr(M,lambda,rho,DEBUG);
%M_hat = M*Z + E;
L = M*Z;
S = E; %S = M_hat - L;