% GreGoDec (Zhou and Tao 2013)
% process_video('RPCA', 'GreGoDec', 'dataset/cctv.avi', 'output/demo_GreGoDec.avi');
rank = 1;
tau = 7;
power = 5;
tol = 1e-3;
k = 1;
L = GreGoDec(M,rank,tau,tol,power,k);
S = M - L;
