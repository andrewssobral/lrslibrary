% Semi-Soft GoDec (Zhou and Tao 2011)
% process_video('RPCA', 'SSGoDec', 'dataset/cctv.avi', 'output/demo_SSGoDec.avi');
rank = 1;
tau = 8;
power = 0;
L = SSGoDec(M,rank,tau,power);
S = M - L;
