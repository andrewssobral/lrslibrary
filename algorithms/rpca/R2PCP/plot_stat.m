%% print/visualize statistics of the solution

line_style='k-';
font_size=16; line_width=1;
% xx=0:k; xx_name='iteration';
xx=stat.cpu; xx_name='CPU time';

figure(102)


subplot(2,3,1)
% plot(0:k,stat.ener,line_style,'LineWidth',line_width),hold on
semilogy(xx,stat.ener,line_style,'LineWidth',line_width),hold on
xlabel(xx_name,'FontSize',font_size)
ylabel('objective value','FontSize',font_size)
set(gca,'FontSize',font_size,'LineWidth',line_width)
% title('objective value')

subplot(2,3,4)
semilogy(xx,stat.f_a,line_style,'LineWidth',line_width),hold on
xlabel(xx_name,'FontSize',font_size)
ylabel('residual norm','FontSize',font_size)
set(gca,'FontSize',font_size,'LineWidth',line_width)
% title('residual norm')

% subplot(2,3,3)
% semilogy(xx,stat.f_b,line_style,'LineWidth',line_width),hold on
% xlabel(xx_name,'FontSize',font_size)
% ylabel('||grad_Bf||','FontSize',font_size)
% set(gca,'FontSize',font_size,'LineWidth',line_width)
% % title('residual norm')

subplot(2,3,2)
semilogy(xx,stat.err_a,line_style,'LineWidth',line_width),hold on
% semilogy(0:k,stat.diff_a,line_style,'LineWidth',line_width),hold on
xlabel(xx_name,'FontSize',font_size)
ylabel('relative error for A','FontSize',font_size)
set(gca,'FontSize',font_size,'LineWidth',line_width)
% title('relative error')

subplot(2,3,3)
semilogy(xx,stat.err_b,line_style,'LineWidth',line_width),hold on
% semilogy(0:k,stat.diff_b,line_style,'LineWidth',line_width),hold on
xlabel(xx_name,'FontSize',font_size)
ylabel('relative error for B','FontSize',font_size)
set(gca,'FontSize',font_size,'LineWidth',line_width)
% title('relative error')

% subplot(2,3,4)
% semilogy(xx(2:end),stat.diff_a,line_style,'LineWidth',line_width),hold on
% xlabel(xx_name,'FontSize',font_size)
% ylabel('||A^k-A^*||','FontSize',font_size)
% set(gca,'FontSize',font_size,'LineWidth',line_width)
% % title('relative difference')
% 
% subplot(2,3,5)
% semilogy(xx(2:end),stat.diff_b,line_style,'LineWidth',line_width),hold on
% xlabel(xx_name,'FontSize',font_size)
% ylabel('||B^k-B^*||','FontSize',font_size)
% set(gca,'FontSize',font_size,'LineWidth',line_width)
% % title('relative difference')

% subplot(2,3,5)
% plot(xx(2:end),stat.ang1,line_style,'LineWidth',line_width),hold on
% xlabel(xx_name,'FontSize',font_size)
% ylabel('\delta_A','FontSize',font_size)
% set(gca,'FontSize',font_size,'LineWidth',line_width)
% % axis([1 k 0 1])
% 
% subplot(2,3,6)
% plot(xx(2:end),stat.ang2,line_style,'LineWidth',line_width),hold on
% xlabel(xx_name,'FontSize',font_size)
% ylabel('\delta_b','FontSize',font_size)
% set(gca,'FontSize',font_size,'LineWidth',line_width)
% % axis([1 k 0 1])

subplot(2,3,5)
plot(xx,stat.n3/n1,line_style,'LineWidth',line_width),hold on
xlabel(xx_name,'FontSize',font_size)
ylabel('rank(A^k)/m','FontSize',font_size)
set(gca,'FontSize',font_size,'LineWidth',line_width)

subplot(2,3,6)
plot(xx,stat.n4/(n1*n2),line_style,'LineWidth',line_width),hold on
xlabel(xx_name,'FontSize',font_size)
ylabel('||B^k||_0/(mn)','FontSize',font_size)
set(gca,'FontSize',font_size,'LineWidth',line_width)


