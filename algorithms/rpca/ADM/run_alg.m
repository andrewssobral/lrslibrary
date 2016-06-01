% ADM / LRSD (Yuan and Yang 2009)
% process_video('RPCA', 'ADM', 'dataset/demo.avi', 'output/demo_ADM.avi');
% TODO: ----> Works only on win32 (mexsvd.mexw32)
t = 0.01;
opts = [];
opts.beta = .25/mean(abs(M(:))); % 0.10;
opts.tol = 1e-6;
opts.maxit = 100; %1000
opts.print = 1;
out = ADM(M, t/(1-t), opts);
L = out.LowRank;
S = out.Sparse;