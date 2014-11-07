function results = horpca_lambda_search( dataname )

runRPCA = false;
results1 = [];  results2 = [];  results3=  [];  results4 = [];  results5 = [];  
results6 = [];  results7 = [];  results7b = [];
lambdaS = 1;

% no missing data
IsTC = false;

rr = [3,2,1];    %[ 1/2, 1/1.8, 1/1.6, 1/1.4, 1/1.2, 1 ];
K = length(rr);
results1 = cell( 1, K );
for i = 1:K
    results1{i} = test_trpca( dataname, 1, rr(i), lambdaS, IsTC, false );
end
disp( 'results1 done.' );

rr = [1,1/2,1/3];    %[1/10, 1/8, 1/6, 1/4, 1/2, 1];
K = length(rr);
results4 = cell( 1, K );
for i = 1:K
    results4{i} = test_trpca( dataname, 4, rr(i), lambdaS, IsTC, false );
end
disp( 'results4 done.' );

rr = [8,6,4,2,1];    %[1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3];
K = length(rr);
results6 = cell( 1, K );
for i = 1:K
    results6{i} = test_trpca( dataname, 6, rr(i), lambdaS, IsTC, false );
end
disp( 'results6 done.' );

rr = [1,1/2,1/3];    %[1/10, 1/8, 1/6, 1/4, 1/2, 1];
K = length(rr);
results7 = cell( 1, K );
for i = 1:K
    results7{i} = test_trpca( dataname, 7, rr(i), 10, IsTC, false );
end
disp( 'results7a done.' );

% with missing data
IsTC = true;

rr = [3,2,1];%[1/1.8, 1/1.6, 1/1.4, 1/1.2, 1, 1/0.8, 1/0.6];
K = length(rr);
results2 = cell( 1, K );
for i = 1:K
    results2{i} = test_trpca( dataname, 2, rr(i), lambdaS, IsTC, false );
end
disp( 'results2 done.' );

rr = [1,1/2,1/3];%[1/10, 1/8, 1/6, 1/4, 1/2, 1];
K = length(rr);
results5 = cell( 1, K );
for i = 1:K
    results5{i} = test_trpca( dataname, 5, rr(i), lambdaS, IsTC, false );
end
disp( 'results5 done.' );

rr = [1,1/2,1/3];%[1/10, 1/8, 1/6, 1/4, 1/2, 1];
K = length(rr);
results7b = cell( 1, K );
for i = 1:K
    results7b{i} = test_trpca( dataname, 7, rr(i), 10, IsTC, false );
end
disp( 'results7b done.' );

results = {};
results{1} = results1;
results{2} = results2;
results{3} = results3;
results{4} = results4;
results{5} = results5;
results{6} = results6;
results{7} = results7;
results{8} = results7b;
end