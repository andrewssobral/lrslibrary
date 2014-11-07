% Demo program for NeNMF and its extensions

%############################
% Test NeNMF
fprintf('Test NeNMF on both synthetic and real datasets.\n');

% synthetic dense dataset
disp('press any key to continue');  pause;
V=max(0,randn(500,100));
reducedDim=50;
fprintf('\nNeNMF factorizing (dense) sub-Gaussian matrix (%dx%d-D) with reduced dimensionality = %d ...\n',size(V,1),size(V,2),reducedDim);
[W,H,it,ela,HIS]=NeNMF(V,reducedDim);
fprintf('Suceeds in %d iterations and %.3f seconds.\n',it,ela);

% synthetic sparse dataset
disp('press any key to continue');  pause;
V=sparse(max(0,randn(500,100)));
reducedDim=50;
fprintf('\nNeNMF factorizing (sparse) sub-Gaussian matrix (%dx%d-D) with reduced dimensionality = %d ...\n',size(V,1),size(V,2),reducedDim);
[W,H,it,ela,HIS]=NeNMF(V,reducedDim);
fprintf('Suceeds in %d iterations and %.3f seconds.\n',it,ela);

% face image dataset
disp('press any key to continue');  pause;
X=load('../Dataset/ORL.mat');
V=X.fea;
V=V/max(V(:));
reducedDim=50;
fprintf('\nNeNMF factorizing (dense) ''ORL'' face image dataset (%dx%d-D) with reduced dimensionality = %d ...\n',size(V,1),size(V,2),reducedDim);
[W,H,it,ela,HIS]=NeNMF(V,reducedDim);
fprintf('Suceeds in %d iterations and %.3f seconds.\n',it,ela);

% document corpus dataset
disp('press any key to continue');  pause;
X=load('../Dataset/TDT2.mat');
V=X.fea;
V=V/max(V(:));
reducedDim=10;
fprintf('\nNeNMF factorizing (sparse) ''TDT2'' corpus dataset (%dx%d-D) with reduced dimensionality = %d ...\n',size(V,1),size(V,2),reducedDim);
[W,H,it,ela,HIS]=NeNMF(V,reducedDim);
fprintf('Suceeds in %d iterations and %.3f seconds.\n\n',it,ela);

%############################
% Test NeNMF Extensions
fprintf('Test NeNMF extensions on both synthetic or real datasets.\n');

% NeNMF-L1 on synthetic dense dataset
disp('press any key to continue');  pause;
V=max(0,randn(500,100));
reducedDim=50;
beta=1;
fprintf('\nNeNMF-L1 factorizing (dense) sub-Gaussian matrix (%dx%d-D) with reduced dimensionality = %d and beta = %.3f...\n',size(V,1),size(V,2),reducedDim,beta);
[W,H,it,ela,HIS]=NeNMF(V,reducedDim,'TYPE','L1R');
fprintf('Suceeds in %d iterations and %.3f seconds.\n',it,ela);

% NeNMF-L2 on synthetic dense dataset
disp('press any key to continue');  pause;
V=max(0,randn(500,100));
reducedDim=50;
beta=1;
fprintf('\nNeNMF-L2 factorizing (dense) sub-Gaussian matrix (%dx%d-D) with reduced dimensionality = %d and beta = %.3f...\n',size(V,1),size(V,2),reducedDim,beta);
[W,H,it,ela,HIS]=NeNMF(V,reducedDim,'TYPE','L2R');
fprintf('Suceeds in %d iterations and %.3f seconds.\n',it,ela);

% NeNMF-MR on real dataset
disp('press any key to continue');  pause;
X=load('../Dataset/ORL.mat');
V=X.fea;    S=constructW(V);
V=V'/max(V(:));
reducedDim=50;
beta=1;
fprintf('\nNeNMF-MR factorizing (dense) ''ORL'' face image dataset (%dx%d-D) with reduced dimensionality = %d and beta = %.3f...\n',size(V,1),size(V,2),reducedDim,beta);
[W,H,it,ela,HIS]=NeNMF(V,reducedDim,'TYPE','MR','S_MTX',S);
fprintf('Suceeds in %d iterations and %.3f seconds.\n',it,ela);