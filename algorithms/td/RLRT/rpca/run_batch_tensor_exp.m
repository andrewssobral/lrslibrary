%function [Err, Cpu, Iter] = run_batch_tensor_exp( dataname )
function [Err, Cpu, Iter] = run_batch_tensor_exp(data)

%load( ['..\..\data\', dataname] );

marker = get_markers();
alg_names = get_algs();

Nrep = 3;   %data.Nrep;
Nobs = length( data.obs );
% algs = [ 2, 5, 7, 8, 6, 3 ];
algs = [2,3];
Nalgs = length(algs);
sz = size( data.X );
N = ndims(data.X);
IsTC = true;
verbose = false;
Err = {};
Cpu = {};
Iter = {};
count = 1;

% vary % observations
Imag = 1;
% 10% noise, mag 1
Inoise = 3;
% rRatios = [ 1.6, 1/4, 1/4, 1, 2, 2 ];
% lambdaS = [ 1, 1, -1, -1, 1, 1 ];
rRatios = [ 1.6, 2 ];
lambdaS = [ 1, 1 ];
error = zeros( Nobs, Nrep, Nalgs );
cpu = zeros( Nobs, Nrep, Nalgs );
iter = zeros( Nobs, Nrep, Nalgs );
for i = 1:Nobs
    if i > 16
%         rRatios = [ 1.5, 1/4, 1/4, 1, 2, 2 ];
%         lambdaS = [ 1, 1, -1, -1, 1, 1 ];
        rRatios = [ 1.5, 2 ];
        lambdaS = [ 1, 1 ];
    end
    for k = 1:Nrep
        dtemp = gen_syn_data( data, i, Inoise, Imag, k );
        for j = 1:Nalgs
            results = test_trpca( dtemp, algs(j), rRatios(j), lambdaS(j), IsTC, verbose, 1 );
            error( i, k, j ) = results.rel_err1;
            cpu( i, k, j ) = results.cpu;
            iter( i, k, j ) = results.iter;
        end
    end
    terr = tensor(mean(error,2));
    double(terr(i,1,:))
end
avg_err = mean(error,2);
Err{count} = avg_err;
avg_cpu = mean(cpu, 2);
Cpu{count} = avg_cpu;
avg_iter = mean(iter,2);
Iter{count} = avg_iter;
figure(count);
subplot( 1, 2, 1 );
for j = 1:Nalgs
    plot( data.obsPct, avg_err(:,1,j)', marker{algs(j)}, 'LineWidth', 2.0 ); hold on;
end
title( ['Error v.s. Observations (noise=', num2str(data.noisePct(Inoise)), ', magnitude=', num2str(data.mag(Imag)), ')'] );
legend( alg_names(algs) );
xlabel( 'Fraction of observations' );
ylabel( 'Relative error' );
hold off;
subplot( 1, 2, 2 );
for j = 1:Nalgs
    plot( data.obsPct, avg_iter(:,1,j)', marker{algs(j)}, 'LineWidth', 2.0 ); hold on;
end
title( ['No. Iters v.s. Observations (noise=', num2str(data.noisePct(Inoise)), ', magnitude=', num2str(data.mag(Imag)), ')'] );
legend( alg_names(algs) );
xlabel( 'Fraction of observations' );
ylabel( 'No. iters' );
hold off;
count = count + 1;



% vary % noise
% 100% obs, mag 1
% algs = [ 1, 4, 7, 8, 6, 3 ];
algs = [1,3];
IsTC = false;
Nnoise = 9;
rRatios = [ 1, 1, 1, 1, 1, 1/1.5, 1/1.5, 1/1.5, 1/1.7; ...
%             1/3, 1/3.5, 1/4, 1/4.5, 1/5, 1/5.5, 1/6, 1/6.5, 1/7; ...
%             1/3, 1/3.5, 1/4, 1/4.5, 1/5, 1/5.5, 1/6, 1/6.5, 1/7; ...
%             1, 1, 1, 1, 1, 1/1.5, 1/1.5, 1/1.5, 1/1.7; ...
%             2, 2, 2, 1.5, 1.5, 1.5, 1, 1, 1; ...
            2, 2, 2, 1.5, 1.5, 1.5, 1, 1, 1];
lambdaS = [ 1, 1, 1, 1, 1, 1, 1, 1, 1; ...
%             1, 1, 1, 1, 1, 1, 1, 1, 1; ...
%             -1, -1, -1, -1, 1/4, 1/4, 1/4, 1/4, 1/4; ...
%             -1, -1, -1, -1, 1/8, 1/8, 1/8, 1/8, 1/8; ...
%             1, 1, 1, 1, 1, 1, 1, 1, 1; ...
            1, 1, 1, 1, 1, 1, 1, 1, 1;];

Iobs = 20;
Imag = 1;
error = zeros( Nnoise, Nrep, Nalgs );
cpu = zeros( Nnoise, Nrep, Nalgs );
iter = zeros( Nnoise, Nrep, Nalgs );
for i = 1:Nnoise
    for k = 1:Nrep
        dtemp = gen_syn_data( data, Iobs, i, Imag, k );     % hardcode k=2 !!!
        for j = 1:Nalgs
            results = test_trpca( dtemp, algs(j), rRatios(j,i), lambdaS(j,i), IsTC, verbose, 1 );
            error( i, k, j ) = results.rel_err1;
            cpu( i, k, j ) = results.cpu;
            iter( i, k, j ) = results.iter;
        end
    end
    terr = tensor(mean(error,2));
    double(terr(i,1,:))
end
avg_err = mean(error,2);
Err{count} = avg_err;
avg_cpu = mean(cpu, 2);
Cpu{count} = avg_cpu;
avg_iter = mean(iter,2);
Iter{count} = avg_iter;
figure(count);
subplot( 1, 2, 1 );
for j = 1:Nalgs
    plot( data.noisePct(1:Nnoise), avg_err(:,1,j)', marker{algs(j)}, 'LineWidth', 2.0 ); hold on;
end
title( ['Error v.s. Percentage Corruptions (obs=', num2str(data.obsPct(Iobs)), ', noise magnitude=', num2str(data.mag(Imag)), ')'] );
legend( alg_names(algs) );
xlabel( 'Percentage corruptions' );
ylabel( 'Relative error' );
hold off;
subplot( 1, 2, 2 );
for j = 1:Nalgs
    plot( data.noisePct(1:Nnoise), avg_iter(:,1,j)', marker{algs(j)}, 'LineWidth', 2.0 ); hold on;
end
title( ['No. Iters v.s. Percentage Corruptions (obs=', num2str(data.obsPct(Iobs)), ', noise magnitude=', num2str(data.mag(Imag)), ')'] );
legend( alg_names(algs) );
xlabel( 'Percentage corruptions' );
ylabel( 'No. iters' );
hold off;
count = count + 1;




end
