function [B, optimal_X, optimal_S, delta]=create_data_noisy_L2(n1,n2,p,r,stdev)

    % random +/- 1 signal
    optimal_S = zeros(n1,n2);
    optimal_X = randn(n1,r)*randn(r,n2);
    noise = stdev*randn(n1,n2);
    %delta = norm(noise(:),2);
    delta = stdev*sqrt(n1*n2);
    OMEGA = randperm(n1*n2);
    OMEGA = OMEGA(1:p);
    optimal_S(OMEGA) = -100+200*rand(1,p);
    B = optimal_X+optimal_S+noise;