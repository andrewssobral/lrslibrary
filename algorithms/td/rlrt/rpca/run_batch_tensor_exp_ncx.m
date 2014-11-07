function [Err, Cpu, Iter] = run_batch_tensor_exp_ncx( dataname )

load( ['..\..\data\', dataname] );

marker = get_markers();
alg_names = get_algs();

Nrep = 3;   %data.Nrep;
Nobs = length( data.obs );
algs = [ 11, 11, 11 ];
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

% 25% noise, mag 1  
algs = [ 11, 11, 11 ];
ks = [ 5 5 5; 6 6 6; 8 8 8 ];
Inoise = 6;
error = zeros( Nobs, Nrep, Nalgs );
cpu = zeros( Nobs, Nrep, Nalgs );
iter = zeros( Nobs, Nrep, Nalgs );
for i = 1:Nobs
    for k = 1:Nrep
        dtemp = gen_syn_data( data, i, Inoise, Imag, k );
        parfor j = 1:Nalgs
%             results = test_trpca( dtemp, algs(j), rRatios(j)/(base0+i*inc), lambdaS(j), IsTC, verbose );
            results = test_trpca( dtemp, algs(j), 1, 1, IsTC, verbose, [], 1, ks(j,:) );
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
% save syn-obs-025-noise-ncx Err Cpu Iter;

figure(count);
subplot( 1, 2, 1 );
for j = 1:Nalgs
    plot( data.obsPct, avg_err(:,1,j)', marker{10+j}, 'LineWidth', 2.0 ); hold on;
end
title( ['Error v.s. Observations (noise=', num2str(data.noisePct(Inoise)), ', magnitude=', num2str(data.mag(Imag)), ')'] );
legend( '5 5 5', '6 6 6', '8 8 8' );
xlabel( 'Fraction of observations' );
ylabel( 'Relative error' );
hold off;
subplot( 1, 2, 2 );
for j = 1:Nalgs
    plot( data.obsPct, avg_iter(:,1,j)', marker{10+j}, 'LineWidth', 2.0 ); hold on;
end
title( ['No. Iters v.s. Observations (noise=', num2str(data.noisePct(Inoise)), ', magnitude=', num2str(data.mag(Imag)), ')'] );
legend( '5 5 5', '6 6 6', '8 8 8' );
xlabel( 'Fraction of observations' );
ylabel( 'No. iters' );
hold off;
count = count + 1;


% vary % noise
% 100% obs, mag 1
algs = [ 10, 10, 10 ];
ks = [ 5 5 5; 6 6 6; 8 8 8 ];
IsTC = false;
Nnoise = 10;
Iobs = 20;
Imag = 1;
error = zeros( Nnoise, Nrep, Nalgs );
cpu = zeros( Nnoise, Nrep, Nalgs );
iter = zeros( Nnoise, Nrep, Nalgs );
for i = 1:Nnoise
    for k = 1:Nrep
        dtemp = gen_syn_data( data, Iobs, i, Imag, k );     % hardcode k=2 !!!
        parfor j = 1:Nalgs
            results = test_trpca( dtemp, algs(j), 1, 1, IsTC, verbose, [], 1, ks(j,:) );
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
% save syn-noise-ncx Err Cpu Iter;

figure(count);
subplot( 1, 2, 1 );
for j = 1:Nalgs
    plot( data.noisePct(1:Nnoise), avg_err(:,1,j)', marker{10+j}, 'LineWidth', 2.0 ); hold on;
end
title( ['Error v.s. Percentage Corruptions (obs=', num2str(data.obsPct(Iobs)), ', noise magnitude=', num2str(data.mag(Imag)), ')'] );
legend( '5 5 5', '6 6 6', '8 8 8' );
xlabel( 'Percentage corruptions' );
ylabel( 'Relative error' );
hold off;
subplot( 1, 2, 2 );
for j = 1:Nalgs
    plot( data.noisePct(1:Nnoise), avg_iter(:,1,j)', marker{10+j}, 'LineWidth', 2.0 ); hold on;
end
title( ['No. Iters v.s. Percentage Corruptions (obs=', num2str(data.obsPct(Iobs)), ', noise magnitude=', num2str(data.mag(Imag)), ')'] );
legend( '5 5 5', '6 6 6', '8 8 8' );
xlabel( 'Percentage corruptions' );
ylabel( 'No. iters' );
hold off;
count = count + 1;

% % partially low-rank data, 10% noise, mag = 1
% % vary % obs
% load ..\..\data\syn-plr
% algs = [ 11, 11, 11 ];
% ks = [ 3 3 20 10; 4 5 20 10; 4 5 17 8 ];
% Nalgs = length(algs);
% Imag = 1;
% Inoise = 3;
% error = zeros( Nobs, Nrep, Nalgs );
% cpu = zeros( Nobs, Nrep, Nalgs );
% iter = zeros( Nobs, Nrep, Nalgs );
% for i = 1:Nobs
%     for k = 1:Nrep
%         dtemp = gen_syn_data( data, i, Inoise, Imag, k );
%         parfor j = 1:Nalgs
%             results = test_trpca( dtemp, algs(j), 1, 1, IsTC, verbose, [], 1, ks(j,:) );
%             error( i, k, j ) = results.rel_err1;
%             cpu( i, k, j ) = results.cpu;
%             iter( i, k, j ) = results.iter;
%         end
%     end
%     terr = tensor(mean(error,2));
%     double(terr(i,1,:))
% end
% avg_err = mean(error,2);
% Err{count} = avg_err;
% avg_cpu = mean(cpu, 2);
% Cpu{count} = avg_cpu;
% avg_iter = mean(iter,2);
% Iter{count} = avg_iter;
% save syn-ncx-plr-avg Err Cpu Iter;
% 
% figure(count);
% subplot( 1, 2, 1 );
% for j = 1:Nalgs
%     plot( data.obsPct, avg_err(:,1,j)', marker{10+j}, 'LineWidth', 2.0 ); hold on;
% end
% title( ['Error v.s. Observations (noise=', num2str(data.noisePct(Inoise)), ', magnitude=', num2str(data.mag(Imag)), ')'] );
% legend( '3 3 20 10', '4 5 20 10', '4 5 17 8' );
% xlabel( 'Fraction of observations' );
% ylabel( 'Relative error' );
% hold off;
% subplot( 1, 2, 2 );
% for j = 1:Nalgs
%     plot( data.obsPct, avg_iter(:,1,j)', marker{10+j}, 'LineWidth', 2.0 ); hold on;
% end
% title( ['No. Iters v.s. Observations (noise=', num2str(data.noisePct(Inoise)), ', magnitude=', num2str(data.mag(Imag)), ')'] );
% legend( '3 3 20 10', '4 5 20 10', '4 5 17 8' );
% xlabel( 'Fraction of observations' );
% ylabel( 'No. iters' );
% hold off;
% count = count + 1;



end
