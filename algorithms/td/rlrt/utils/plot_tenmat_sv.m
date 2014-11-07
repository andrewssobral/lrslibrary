function Nsv = plot_tenmat_sv( T, fig_ind )
% T is a tensor
% plots the SV's of unfoldings in each mode
% return Nsv = no. of SV's larger than 0.01 * max SV


if ~exist('fig_ind', 'var') || isempty(fig_ind)
    fig_ind = 10;
end

N = length(size(T));
figure(fig_ind);
scrsz = get(0,'ScreenSize');
set(fig_ind,'Position',[10 scrsz(4)*0.05 scrsz(3)*0.4 0.8*scrsz(4)]);
Nsv = zeros( 1, N );

for mode = 1:N
    s = svd( double(tenmat(T,mode)) );
    sMax = max(s);
    Isv = s > 0.01*sMax;        % before: 0.01!!!
    Nsv(mode) = sum( Isv );
    
    subplot( N, 1, mode );
    plot( s, 'LineWidth', 2 );
    title( ['Singular values for mode ', num2str(mode)] );
%     subplot( 2, N, mode+N );
%     s(~Isv) = 0;
%     plot( s, 'LineWidth', 2 );
end

end