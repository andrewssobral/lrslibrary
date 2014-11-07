%% to normalizes the values of A to values between 0 and 1;

function B=imnorm(A)

min_A=min(A(:));
max_A=max(A(:));

B= (A-min_A)/(max_A-min_A);

