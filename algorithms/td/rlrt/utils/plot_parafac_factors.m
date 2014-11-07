function [P_orig, P] = plot_parafac_factors( results, K, Ko, isSing, fig_ind )

P_orig = parafac_als( results.X0, Ko, struct('maxiters',100) );
P = horpca_parafac( results, K, isSing );
N = length(size(results.X));

if ~exist('fig_ind', 'var') || isempty(fig_ind)
    fig_ind = 10;
end
figure(fig_ind);
scrsz = get(0,'ScreenSize');
set(fig_ind,'Position',[10 scrsz(4)*0.05 scrsz(3)*0.7 0.8*scrsz(4)]);

for i = 1:N
    subplot( N, 2, 1+(i-1)*2 );
    plot( P_orig.U{i}, 'LineWidth', 2 );
    title( ['Original loading factors for mode ', num2str(i)] );
    
    subplot( N, 2, 2+(i-1)*2 );
    plot( P.U{i}, 'LineWidth', 2 );
    title( ['Reconstructed loading factors for mode ', num2str(i)] );
end

end