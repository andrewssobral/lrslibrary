% Markov BRMF (Wang and Yeung 2013)
% process_video('RPCA', 'MBRMF', 'dataset/demo.avi', 'output/demo_MBRMF.avi');
D = normalize(M);
% Set up the (hyper)parameters. See more details in the paper.
% r = 10;
r = round(min(20,sqrt(size(D,2)))/2);
opts.maxIter = 100;
opts.burnin = 50;
opts.invW_0 = 1000 * eye(r * 2);
opts.beta_0 = 2;
opts.nu_0 = r * 2;
opts.a = 1e-4;
opts.b = 1e0; %[1 ~ 10]
% We set the maximum rank to be twice of the ground truth.
opts.r = r * 2;
opts.alpha = 0.5;
[~,~,~,p] = MBRMF(D, opts);
L = p;
S = (D - p);