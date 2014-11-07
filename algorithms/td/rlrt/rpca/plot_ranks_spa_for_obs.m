function stats = plot_ranks_spa_for_obs( fig_ind, saveResults )
% plot ranks vs corruption % for a few obs

load ..\..\data\syn-lr-tnames   % load tnames

marker = { 'r*', 'bo', 'k+', 'g#' };
Ntens = length(tnames);
% Itens = [ 1,2,3 ];      
% Ntens = length(Itens);
Iobs = [ 20, 16, 12 ];
Nobs = length( Iobs );
Nnoise = 9;
Nrep = 1;
Imag = 1;
alg = 2;
rRatio = 1;
lambdaS = 1;
IsTC = true;
verbose = false;
err_bar = 0.01;

stats.ranks = cell( 1, Ntens );
stats.tranks = zeros( 2, Ntens );
stats.spa = zeros( Nobs, Ntens );
stats.exp = 'rank-spa-for-obs';

for t = 1:Ntens
    tname = tnames{t};
    load( ['..\..\data\',tname] );
    maxNoise = Nnoise;
    stats.ranks{t} = data.rank;
    [ stats.tranks(1,t), stats.tranks(2,t) ] = compute_tranks( size(data.X), data.rank );
    
    for i = 1:Nobs
        obsInd = Iobs(i);
        for j = maxNoise:-1:1
            err = zeros( 1, Nrep );
            for k = 1:Nrep
                dtemp = gen_syn_data( data, obsInd, j, Imag, k );
                results = test_trpca( dtemp, alg, rRatio, lambdaS, IsTC, verbose );
                err(k) = results.rel_err1;
            end
            avg_err = mean(err);
            fprintf( '%s,  obs: %1.2f,  spa: %1.2f,  err: %1.2e\n', tname, data.obsPct(obsInd), data.noisePct(j), avg_err );
            if avg_err <= err_bar
                % record spa
                stats.spa( i, t ) = data.noisePct(j);
                break;  % found the largest corruption % tolerated
            end
        end
        
        maxNoise = j;
        if j == 1 && stats.spa(i,t) == 0    % no noise tolerated from now on
            break;
        end
    end
    
%     plot( stats.obs, stats.spa(t,:), marker{t}, 'LineWidth', 2 );
%     hold on;
end

figure( fig_ind );
for i = 1:Nobs
    plot( stats.tranks(1,:), stats.spa(i,:), marker{i}, 'LineWidth', 2 );
    hold on;
end
title( 'rank sums v.s. max corruption % for various obs %' );
xlabel( 'rank sum' );
ylabel( 'max corruption %' );
legend( numarray2strarray(data.obsPct(Iobs)) );
hold off;

figure( fig_ind*2 );
for i = 1:Nobs
    plot( stats.tranks(2,:), stats.spa(i,:), marker{i}, 'LineWidth', 2 );
    hold on;
end
title( 'normalized ranks v.s. max corruption % for various obs %' );
xlabel( 'normalized rank' );
ylabel( 'max corruption %' );
legend( numarray2strarray(data.obsPct(Iobs)) );
hold off;

if exist( 'saveResults', 'var' ) && ~isempty(saveResults) && saveResults
    save ..\..\data\syn-ranks-spa stats
end

end