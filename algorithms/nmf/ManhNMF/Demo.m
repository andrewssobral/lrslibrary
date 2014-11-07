% Demo for ManhNMF algorithms including RRI and AGD.
% Please include 'ShowImage.m' together with this demo.

% Face Image Illumination Modeling
load ../Dataset/YALEB.mat
X=fea(gnd==1,:)'/255;
[m,n]=size(X);
r=2;

W0=rand(r,m)/r;
H0=rand(r,n)/r;
W0=W0.*(H0*X')./(H0*H0'*W0);
H0=H0.*(W0*X)./(W0*W0'*H0);

fprintf('Illumination Modeling by ManhNMF on Yale B dataset ...\n');
[W,H]=ManhNMF(X,r,'w_init',W0,'h_init',H0,'lam_init',.1,'tol_innr',1e-3,'vb_outr',2);
Z=W'*H;
S=X-Z;

figure(1);
subplot(1,3,1);
ShowImage(X,32,32,8,8,0);
subplot(1,3,2);
Z=Z./(ones(m,1)*max(Z));
ShowImage(Z,32,32,8,8,0);
subplot(1,3,3);
ShowImage(S,32,32,8,8,0);

fprintf('Illumination modeling successes.\n\n');

% Face Image Representation
load ../Dataset/PIE.mat;
[tr,vd,ts]=RandPartDB(gnd,2);
X=fea(tr,:)'/255;
[m,n]=size(X);
r=64;

W0=rand(r,m)/r;
H0=rand(r,n)/r;
W0=W0.*(H0*X')./(H0*H0'*W0);
H0=H0.*(W0*X)./(W0*W0'*H0);

fprintf('AGD algorithm for ManhNMF on PIE dataset ...\n');
[Wagd,Hagd,~,~,HISagd]=ManhNMF(X,r,'w_init',W0,'h_init',H0,...
    'alg_type','agd','vb_innr',0,'vb_outr',2);
Hagd=Hagd.*(sum(Wagd,2)*ones(1,n));
Wagd=Wagd./(sum(Wagd,2)*ones(1,m));

fprintf('RRI algorithm for ManhNMF on PIE dataset ...\n');
[Wrri,Hrri,~,~,HISrri]=ManhNMF(X,r,'w_init',W0,'h_init',H0,...
    'alg_type','rri','vb_innr',0,'vb_outr',2);
Hrri=Hrri.*(sum(Wrri,2)*ones(1,n));
Wrri=Wrri./(sum(Wrri,2)*ones(1,m));

figure(2);
subplot(2,2,1);
semilogy(HISrri.objf,'-ob','LineWidth',2.5,'MarkerSize',8);
hold on;
semilogy(HISagd.objf,'-pr','LineWidth',2.5,'MarkerSize',8);
hold off;
title('(a) objective values vs. iteration number','FontSize',16);
xlabel('iteration numbers','FontSize',16);
ylabel('objective values','FontSize',16);
legend('RRI','AGD');
set(gca,'FontSize',16)

subplot(2,2,2);
semilogy(HISrri.cpus,HISrri.objf,'-ob','LineWidth',2.5,'MarkerSize',8);
hold on;
semilogy(HISagd.cpus,HISagd.objf,'-pr','LineWidth',2.5,'MarkerSize',8);
hold off;
title('(b) objective values vs. CPU seconds','FontSize',16);
xlabel('CPU seconds','FontSize',16);
ylabel('objective values','FontSize',16);
legend('RRI','AGD');
set(gca,'FontSize',16)

subplot(2,2,3);
ShowImage(Wrri',32,32,8,8,0);
title('(c) basis images obtained by RRI for ManhNMF','FontSize',16);
set(gca,'FontSize',16)

subplot(2,2,4);
ShowImage(Wagd',32,32,8,8,0);
title('(d) basis images obtained by AGD for ManhNMF','FontSize',16);
set(gca,'FontSize',16)

fprintf('Face image representation successes.\n\n');