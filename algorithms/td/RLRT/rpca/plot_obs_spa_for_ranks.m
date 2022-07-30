function stats = plot_obs_spa_for_ranks( fig_ind, saveResults )
% plot obs % vs corruption % for a few ranks

load ..\..\data\syn-lr-tnames   % load tnames

marker = { 'r-*', 'b-*', 'k-*', 'g-*' };
Ntens = length(tnames);
Itens = [ 1,2,3 ];      
Ntens = length(Itens);
Nobs = 20;
Nnoise = 9;
Nrep = 1;
Imag = 1;
alg = 2;
rRatio = 1;
lambdaS = 1;
IsTC = true;
verbose = false;
err_bar = 0.01;

stats.obs = [];
stats.spa = zeros( Ntens, Nobs );
stats.exp = 'obs-spa-for-ranks';

figure( fig_ind );

for t = 1:Ntens
    tname = tnames{Itens(t)};
    load( ['..\..\data\',tname] );
    maxNoise = Nnoise;
    
    for i = Nobs:-1:1
        
        for j = maxNoise:-1:1
            err = zeros( 1, Nrep );
            for k = 1:Nrep
                dtemp = gen_syn_data( data, i, j, Imag, k );
                results = test_trpca( dtemp, alg, rRatio, lambdaS, IsTC, verbose );
                err(k) = results.rel_err1;
            end
            avg_err = mean(err);
            fprintf( '%s,  obs: %1.2f,  spa: %1.2f,  err: %1.2e\n', tname, data.obsPct(i), data.noisePct(j), avg_err );
            if avg_err <= err_bar
                % record spa
                stats.spa( t, i ) = data.noisePct(j);
                break;  % found the largest corruption % tolerated
            end
        end
        
        maxNoise = j;
        if j == 1 && stats.spa(t,i) == 0    % no noise tolerated from now on
            break;
        end
    end
    
    stats.obs = data.obsPct(1:Nobs);
    plot( stats.obs, stats.spa(t,:), marker{t}, 'LineWidth', 2 );
    hold on;
end

title( 'obs % v.s. max corruption % for various ranks' );
xlabel( 'obs %' );
ylabel( 'max corruption %' );
legend( tnames(1:Ntens) );
hold off;

if exist( 'saveResults', 'var' ) && ~isempty(saveResults) && saveResults
    save ..\..\data\syn-obs-spa stats
end

end