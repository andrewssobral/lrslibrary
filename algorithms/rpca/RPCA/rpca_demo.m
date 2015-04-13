
% This Software computes the Robust Principal Component Analysis or Robust Singular Value Decomposition
% explained at the references:

% REFERENCES:
% De la Torre, F. and Black, M. J., Robust principal component analysis for computer vision. Int. Conf. on Computer Vision 2001, Vancouver. 

%This software is provided for research purposes only, any sale or use for commercial purposes is strictly prohibited. 

%This is EXPERIMENTAL software, so with high probability it is possible that it has bugs. If you find any bug, please
%send an e-mail to  ftorre@salleURL.edu .
%Use at your own risk. No warranty is implied by this distribution. 

% Copyright 2001. Fernando De la Torre.

% Attention!!!!!
% All the operations are optimized vectorized in order to take profit of the Matlab structure, but to run
% this program without writing in the hard disk requires at least 800MB of RAM memory.

clear all
close all

%load original data
load data_2_days 
% to visualize each of the image data simply write :  imagesc(reshape(Data(:,i),sizeim))
show_2dvideo(Data,120,160);

%Important!!!! if your computer does not have more than 512Mb RAM memory, do 
% Data=Data(:,1:end/2) % it just computes RPCA on the data of 1 day.

%computing standard PCA
mean_ls=mean(Data')';
[u,s,v]=svd(Data-mean_ls*ones(1,size(Data,2)),0);
u=u(:,1:ceil(size(Data,2)/10));s=diag(s);v=v(:,1:ceil(size(Data,2)/10));

%load off-line PCA
%load pca_result_2_days 

Nl=1; % Spatial Neighborhood to compute scale statistics
beta=2.3; % Parameter ....
number_component=8; %number of components which preserve 55% of the energy with standard PCA.

%compute scale statistics
Sigmaf=zeros(sizeim);
cini=(u(:,1:number_component)'*(Data-mean_ls*ones(1,size(Data,2))));
error=Data-mean_ls*ones(1,size(Data,2))-u(:,1:number_component)*cini;
medianat=median(abs(error(:)));
error2=error(:)-medianat;
Sigmaft=sqrt(3)*1.4826*median(abs(error2(:)));

for i=Nl+1:sizeim(1)-Nl-1
      for j=Nl+1:sizeim(2)-Nl-1
         [y,x]=meshgrid(i-Nl:i+Nl,j-Nl:j+Nl);
         ind=sub2ind(sizeim,y(:),x(:));
         errorlo=error(ind,:);
         medianat=median(abs(errorlo(:)));
   		error2=error(ind,:)-medianat;
		   Sigmaf(i,j)=beta*sqrt(3)*1.4826*median(abs(error2(:)));         
      end
end
Sigmaf=Sigmaf(:);
Sigmaf=max(Sigmaf,Sigmaft*ones(size(Sigmaf)));
Sigmai=3*Sigmaf;

clear error error2 medianat Sigmaft errorlo x y i ind index j 

%Initialize RPCA
mean_ini=median(Data')';
basis_ini=u(:,1:number_component)+randn(size(Data,1),number_component);
c_ini=pinv(basis_ini)*(Data-mean_ini*ones(1,size(Data,2)));
%Here we use the principal components as initial guess, 
%but we can add random noise to the basis and run several times the algorithm to get the solution with lowest minimum.


%Compute RPCA
[rob_mean,Bg,cg,info,Sigma]=rob_pca(Data,number_component,300,1,Sigmaf,Sigmai,2,basis_ini,c_ini,mean_ini);   

%Compute outliers
error=Data-rob_mean*ones(1,size(Data,2))-Bg*cg;
Sigmat=(Sigmaf*ones(1,size(Data,2)));
temp=abs(error)<(Sigmat/sqrt(3));
W2=(Sigmat./((Sigmat+error.^2).^2)).*temp;

clear Sigmat temp error


%Add bases until the preserve energy as defined in the paper is bigger than 0.9 (in this case 24 bases)
%to compute weighted principal component analysis, (We use alternated weighted least squares the normalized gradient version is coming soon... )
%The initial guess is the one given by rob_pca, as before
%we can add random noise to it, run several times the algorithm and get the solution with lowest minimum.

bases=24; 
[Bgw,cgw,info,meanw]=weighted_pca(Data,W2,bases,25,2,[Bg 0.0001*randn(size(Bg,1),bases-size(Bg,2))],[cg ; 0.0001*randn(bases-size(Bg,2),size(Data,2))],rob_mean);   

%The images will be reconstructed as   	 rec=meanw*ones(1,size(Data,2))+Bgw*cgw;
%To visualize them just write
%imagesc(reshape(rec(:,1),sizeim)),colormap('gray')
show_2dvideo(rec,120,160);
