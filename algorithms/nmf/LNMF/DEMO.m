% Demo Testing for NMF, LNMF, DNMF, NDLA, DLA, PCA, and FLDA in terms of
% classification
% Written by Naiyang Guan (ny.guan@gmail.com)

load ../Dataset/UMIST.mat;
occFea=RandRandOcc(fea,40,40,20);
[tr,vd,ts]=RandPartDB(gnd,7);
reducedDim=80;

%%% 1. Training NMF, LNMF, DNMF, NDLA on UMIST dataset
fprintf('-----------------------------------------------------\n');
fprintf('Training NMF, PCA, FLDA, and DLA algorithms:\n');
%-----------------------------------------------------
% NDLA
[W_ndla,H_ndla,ela_ndla]=NDLA(fea(tr,:)/255,gnd(tr),[2,15],reducedDim,63.1,0);         % Training
fprintf('NDLA completed in %f seconds.\n',ela_ndla);

%-----------------------------------------------------
% NMF
[W_nmf,H_nmf,ela_nmf]=NMF(fea(tr,:)'/255,reducedDim,0);
fprintf('NMF completed in %f seconds.\n',ela_nmf);

%-----------------------------------------------------
% LNMF
option.verbose=0;
[W_lnmf,H_lnmf,out]=LNMF(fea(tr,:)'/255,reducedDim,option);
ela_lnmf=out.ela;
fprintf('LNMF completed in %f seconds.\n',ela_lnmf);

%-----------------------------------------------------
% DNMF
option.verbose=0;
option.gamma=10;
option.delta=0.01;
[W_dnmf,H_dnmf,out]=DNMF(fea(tr,:)'/255,gnd(tr),reducedDim,option);
ela_dnmf=out.ela;
fprintf('DNMF completed in %f seconds.\n',ela_dnmf);

%-----------------------------------------------------
% PCA
[W_pca,egv_pca,ela_pca]=PCA(fea(tr,:));
fprintf('PCA completed in %f seconds.\n',ela_pca);

%-----------------------------------------------------
% FLDA
opt.Fisherface = 1;
[W_lda,egv_lda,out]=LDA(gnd(tr),opt,fea(tr,:));
ela_lda=out.timeAll;
fprintf('FLDA completed in %f seconds.\n',ela_lda);

%-----------------------------------------------------
% DLA
tr_data=fea(tr,:)*W_pca(:,1:min(reducedDim,size(W_pca,2)));
tic;
W_dla=DLA(tr_data,gnd(tr),2,2,.5);
ela_dla=toc;
fprintf('DLA completed in %f seconds.\n',ela_dla);

%%% 2. Testing on clear dataset
fprintf('-----------------------------------------------------\n');
fprintf('Testing on clear dataset:\n');
%-----------------------------------------------------
% NDLA
gnd_ts=GetLabel(pinv(W_ndla)*fea(tr,:)',pinv(W_ndla)*fea(ts,:)',gnd(tr));     % Testing
accuracy=sum(gnd_ts(:,1)==gnd(ts))/length(ts);
fprintf('NDLA(%d) accuracy=%f.\n',reducedDim,accuracy);

%-----------------------------------------------------
% NMF
gnd_ts=GetLabel(pinv(W_nmf)*fea(tr,:)',pinv(W_nmf)*fea(ts,:)',gnd(tr));
accuracy=sum(gnd_ts(:,1)==gnd(ts))/length(ts);
fprintf('NMF(%d) accuracy=%f.\n',reducedDim,accuracy);

%-----------------------------------------------------
% LNMF
gnd_ts=GetLabel(pinv(W_lnmf)*fea(tr,:)',pinv(W_lnmf)*fea(ts,:)',gnd(tr));
accuracy=sum(gnd_ts(:,1)==gnd(ts))/length(ts);
fprintf('LNMF(%d) accuracy=%f.\n',reducedDim,accuracy);

%-----------------------------------------------------
% DNMF
gnd_ts=GetLabel(pinv(W_dnmf)*fea(tr,:)',pinv(W_dnmf)*fea(ts,:)',gnd(tr));
accuracy=sum(gnd_ts(:,1)==gnd(ts))/length(ts);
fprintf('DNMF(%d) accuracy=%f.\n',reducedDim,accuracy);

%-----------------------------------------------------
% PCA
gnd_ts=GetLabel(W_pca(:,1:min(reducedDim,size(W_pca,2)))'*fea(tr,:)',W_pca(:,1:min(reducedDim,size(W_pca,2)))'*fea(ts,:)',gnd(tr));
accuracy=sum(gnd_ts(:,1)==gnd(ts))/length(ts);
fprintf('PCA(%d) accuracy=%f.\n',reducedDim,accuracy);

%-----------------------------------------------------
% FLDA
gnd_ts=GetLabel(W_lda(:,1:min(reducedDim,size(W_lda,2)))'*fea(tr,:)',W_lda(:,1:min(reducedDim,size(W_lda,2)))'*fea(ts,:)',gnd(tr));
accuracy=sum(gnd_ts(:,1)==gnd(ts))/length(ts);
fprintf('FLDA(%d) accuracy=%f.\n',reducedDim,accuracy);

%-----------------------------------------------------
% DLA
tr_data=fea(tr,:)*W_pca(:,1:min(reducedDim,size(W_pca,2)));
ts_data=fea(ts,:)*W_pca(:,1:min(reducedDim,size(W_pca,2)));
gnd_ts=GetLabel(W_dla(:,1:min(reducedDim,size(W_dla,2)))'*tr_data',W_dla(:,1:min(reducedDim,size(W_dla,2)))'*ts_data',gnd(tr));
accuracy=sum(gnd_ts(:,1)==gnd(ts))/length(ts);
fprintf('DLA(%d) accuracy=%f.\n',reducedDim,accuracy);

%%% 3. Testing on occluded dataset
fprintf('-----------------------------------------------------\n');
fprintf('Testing on occluded dataset:\n');
%-----------------------------------------------------
% NDLA
gnd_ts=GetLabel(pinv(W_ndla)*fea(tr,:)',pinv(W_ndla)*occFea(ts,:)',gnd(tr));  % Testing on occluded set
accuracy=sum(gnd_ts(:,1)==gnd(ts))/length(ts);
fprintf('NDLA(%d,occ) accuracy=%f.\n',reducedDim,accuracy);

%-----------------------------------------------------
% NMF
gnd_ts=GetLabel(pinv(W_nmf)*fea(tr,:)',pinv(W_nmf)*occFea(ts,:)',gnd(tr)); % Testing on occluded set
accuracy=sum(gnd_ts(:,1)==gnd(ts))/length(ts);
fprintf('NMF(%d,occ) accuracy=%f.\n',reducedDim,accuracy);

%-----------------------------------------------------
% LNMF
gnd_ts=GetLabel(pinv(W_lnmf)*fea(tr,:)',pinv(W_lnmf)*occFea(ts,:)',gnd(tr)); % Testing on occluded set
accuracy=sum(gnd_ts(:,1)==gnd(ts))/length(ts);
fprintf('LNMF(%d,occ) accuracy=%f.\n',reducedDim,accuracy);

%-----------------------------------------------------
% DNMF
gnd_ts=GetLabel(pinv(W_dnmf)*fea(tr,:)',pinv(W_dnmf)*occFea(ts,:)',gnd(tr)); % Testing on occluded set
accuracy=sum(gnd_ts(:,1)==gnd(ts))/length(ts);
fprintf('DNMF(%d,occ) accuracy=%f.\n',reducedDim,accuracy);

%-----------------------------------------------------
% PCA
gnd_ts=GetLabel(W_pca(:,1:min(reducedDim,size(W_pca,2)))'*fea(tr,:)',W_pca(:,1:min(reducedDim,size(W_pca,2)))'*occFea(ts,:)',gnd(tr));
accuracy=sum(gnd_ts(:,1)==gnd(ts))/length(ts);
fprintf('PCA(%d,occ) accuracy=%f.\n',reducedDim,accuracy);

%-----------------------------------------------------
% FLDA
gnd_ts=GetLabel(W_lda(:,1:min(reducedDim,size(W_lda,2)))'*fea(tr,:)',W_lda(:,1:min(reducedDim,size(W_lda,2)))'*occFea(ts,:)',gnd(tr));
accuracy=sum(gnd_ts(:,1)==gnd(ts))/length(ts);
fprintf('FLDA(%d,occ) accuracy=%f.\n',reducedDim,accuracy);

%-----------------------------------------------------
% DLA
tr_data=fea(tr,:)*W_pca(:,1:min(reducedDim,size(W_pca,2)));
ts_data=occFea(ts,:)*W_pca(:,1:min(reducedDim,size(W_pca,2)));
gnd_ts=GetLabel(W_dla(:,1:min(reducedDim,size(W_dla,2)))'*tr_data',W_dla(:,1:min(reducedDim,size(W_dla,2)))'*ts_data',gnd(tr));
accuracy=sum(gnd_ts(:,1)==gnd(ts))/length(ts);
fprintf('DLA(%d,occ) accuracy=%f.\n',reducedDim,accuracy);