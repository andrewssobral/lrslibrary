function stats = plot_ranks_obs_for_spa( fig_ind, saveResults )
% plot ranks vs corruption % for a few obs

load ..\..\data\syn-lr-tnames   % load tnames

marker = { 'r*', 'bo', 'k+', 'g#' };
Ntens = length(tnames);
Nobs = 20;
Inoise = [ 1, 3 ];
Nnoise = length(Inoise);
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
stats.obs = ones( Nnoise, Ntens );
stats.exp = 'rank-obs-for-spa';

for t = 1:Ntens
    tname = tnames{t};
    load( ['..\..\data\',tname] );
    minObs = 1;
    stats.ranks{t} = data.rank;
    [ stats.tranks(1,t), stats.tranks(2,t) ] = compute_tranks( size(data.X), data.rank );
    
    for j = 1:Nnoise
        spaInd = Inoise(j);
        for i = minObs:Nobs
            err = zeros( 1, Nrep );
            for k = 1:Nrep
                dtemp = gen_syn_data( data, i, spaInd, Imag, k );
                results = test_trpca( dtemp, alg, rRatio, lambdaS, IsTC, verbose );
                err(k) = results.rel_err1;
            end
            avg_err = mean(err);
            fprintf( '%s,  obs: %1.2f,  spa: %1.2f,  err: %1.2e\n', tname, data.obsPct(i), data.noisePct(spaInd), avg_err );
            if avg_err <= err_bar
                stats.obs( j, t ) = data.obsPct(i);
                break;  % found the largest corruption % tolerated
            end
        end
        
        minObs = i;
        if i == Nobs && stats.obs(j,t) == 1    % full obs needed from now on
            break;
        end
    end
    
end

figure( fig_ind );
for j = 1:Nnoise
    plot( stats.tranks(1,:), stats.obs(j,:), marker{j}, 'LineWidth', 2 );
    hold on;
end
title( 'rank sums v.s. min obs % for various corruptions %' );
xlabel( 'rank sum' );
ylabel( 'min obs %' );
legend( numarray2strarray(data.noisePct(Inoise)) );
hold off;

figure( fig_ind*2 );
for j = 1:Nnoise
    plot( stats.tranks(2,:), stats.obs(j,:), marker{j}, 'LineWidth', 2 );
    hold on;
end
title( 'normalized ranks v.s. min obs % for various corruptions %' );
xlabel( 'normalized rank' );
ylabel( 'min obs %' );
legend( numarray2strarray(data.noisePct(Inoise)) );
hold off;

if exist( 'saveResults', 'var' ) && ~isempty(saveResults) && saveResults
    save ..\..\data\syn-ranks-obs stats
end

end