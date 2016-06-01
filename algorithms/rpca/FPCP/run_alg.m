% Fast PCP (Rodriguez and Wohlberg 2013)
% process_video('RPCA', 'FPCP', 'dataset/demo.avi', 'output/demo_FPCP.avi');
lambda = 1/sqrt(max(size(M))); % default lambda
[L,S] = fastpcp(M,lambda);