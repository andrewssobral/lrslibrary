% TTD | ADMM | Alternating Direction Method of Multipliers (Parikh and Boyd, 2014)
% process_video('TTD', 'ADMM', 'dataset/demo.avi', 'output/demo_ADMM.avi');

results = ADMM(M); % show_2dvideo(M,m,n);
%E = results.Z; % show_2dvideo(E,m,n);
S = results.S; % show_2dvideo(S,m,n);
L = results.L; % show_2dvideo(L,m,n);
%M_hat = L + S + E; % show_2dvideo(M_hat,m,n);