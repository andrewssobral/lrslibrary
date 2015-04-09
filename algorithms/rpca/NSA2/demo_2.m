seed = 601;
randn('state',seed); rand('state',seed);
global D X S
data_type = 0; frame_array=[35,100,125];
DB=50;
tol = 1e-2;

n1 = 500; n2=500;
p = ceil(n1*n2*0.05);
r = ceil(min(n1,n2)*0.05);
pwr = r+p/(n1*n2)*100^2/3;
stdev = sqrt(pwr/10^(DB/10))

[D, optimal_X, optimal_S, deltabar]=create_data_noisy_L2(n1,n2,p,r,stdev);
[X,S]=nsa_v1(D,stdev,tol,optimal_X,optimal_S);

Nsub = max(n1,n2);
comp = randperm(n1*n2);
comp = comp(1:Nsub);
    
figure;
subplot('Position',[0.1,0.77,0.55,0.2]);
plot([1:Nsub]',optimal_S(comp),'r-','LineWidth',0.75);
set(gca,'XLim',[0,Nsub],'YLim',[min(optimal_S(:)),max(optimal_S(:))],'YTick',[min(optimal_S(:)),max(optimal_S(:))],'XTick',[0:1500:Nsub]);
ylabel('S^0_i','FontSize',14);
subplot('Position',[0.1,0.50,0.55,0.2]);
plot([1:Nsub]',S(comp),'b-','LineWidth',0.75);
set(gca,'XLim',[0,Nsub],'YLim',[min(S(:)),max(S(:))],'YTick',[min(S(:)),max(S(:))],'XTick',[0:1500:Nsub]);
ylabel('S^{sol}_i','FontSize',14);
%ind=find(abs(optimal_S(:))>0);
%min_mag = 10*min(abs(optimal_S(ind)));
min_mag = 0.5;
%min_mag = 1e-2;

subplot('Position',[0.1,0.28,0.55,0.15]);
plot([1:Nsub]',abs(S(comp)-optimal_S(comp)),'b-','LineWidth',0.75);
set(gca,'XLim',[0,Nsub],'YLim',[0,min_mag],'YTick',[0,max(abs(S(comp)-optimal_S(comp))),min_mag],'XTick',[0:1500:Nsub]);
ylabel('|S^{sol}_i-S^0_i|','FontSize',14);

Err = D-optimal_X-optimal_S;
subplot('Position',[0.1,0.06,0.55,0.15]);
plot([1:Nsub]',abs(Err(comp)),'b-','LineWidth',0.75);
set(gca,'XLim',[0,Nsub],'YLim',[0,min_mag],'YTick',[0,max(abs(Err(comp))),min_mag],'XTick',[0:1500:Nsub]);
ylabel('|\zeta^0_i|','FontSize',14);
xlabel('i: components','FontSize',14);    

Sigma_opt = svd(optimal_X, 'econ');
Sigma_X = svd(X,'econ');   
Nsub=100;
subplot('Position',[0.75,0.77,0.15,0.2]);
plot([1:Nsub]',Sigma_opt(1:Nsub),'r.','LineWidth',0.75);
set(gca,'XLim',[0,Nsub],'YLim',[0,max(Sigma_opt(:))],'YTick',[0,Sigma_opt(r),max(Sigma_opt(:))],'XTick',[0:1500:Nsub]);
ylabel('\sigma_i(X^0)','FontSize',14);
subplot('Position',[0.75,0.50,0.15,0.2]);
plot([1:Nsub]',Sigma_X(1:Nsub),'b.','LineWidth',0.75);
set(gca,'XLim',[0,Nsub],'YLim',[0,max(Sigma_X(:))],'YTick',[0,Sigma_X(r),max(Sigma_X(:))],'XTick',[0:1500:Nsub]);
ylabel('\sigma_i(X^{sol})','FontSize',14);
subplot('Position',[0.75,0.28,0.15,0.15]);
plot([1:Nsub]',abs(Sigma_X(1:Nsub)-Sigma_opt(1:Nsub)),'b-','LineWidth',0.75);
%set(gca,'XLim',[0,Nsub],'YLim',[0,Sigma_opt(r)],'YTick',[0,Sigma_opt(r)],'XTick',[0:1500:Nsub]);
set(gca,'XLim',[0,Nsub],'YLim',[0,10*min_mag],'YTick',[0,max(abs(Sigma_X(1:Nsub)-Sigma_opt(1:Nsub))),10*min_mag],'XTick',[0:1500:Nsub]);
ylabel('|\sigma_i(X^{sol})-\sigma_i(X^0)|','FontSize',14);
xlabel('i: components','FontSize',14);
clear all
