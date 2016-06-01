% GoDec (Zhou and Tao 2011)
% process_video('RPCA', 'GoDec', 'dataset/cctv.avi', 'output/demo_GoDec.avi');
rank = 1;
card = numel(M); %card = 3.1e+5;
power = 0;
[L,S] = GoDec(M,rank,card,power);
