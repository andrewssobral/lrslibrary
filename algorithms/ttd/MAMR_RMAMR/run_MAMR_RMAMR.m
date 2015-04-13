%{
clc;
lrs_load_conf;
load(fullfile(lrs_conf.lrs_dir,'dataset','trafficdb','traffic_patches.mat'));
V = im2double(imgdb{100});
[M,m,n,p] = convert_video3d_to_2d(V);
%}

D = M;
L = zeros(size(D));
S = zeros(size(D));
E = zeros(size(D));

%CM = ones(size(D));
%CM = zeros(size(D));
[numr,numc] = size(D);
CM = randi([0 1],numr,numc); % simulated confidence map (binary matrix)

Omega = find(CM ~= 0);
[I, J] = ind2sub([numr numc],Omega);  

% Motion-Assisted Matrix Restoration (Ye et al. 2015)
if(strcmp(algorithm_id,'MAMR'))
  lambda = 1; 
  [L,S,iter1] = core_MAMR(D,lambda,I,J);
end

% Robust Motion-Assisted Matrix Restoration (Ye et al. 2015)
if(strcmp(algorithm_id,'RMAMR'))
  lambda = 10; 
  [L,S,E,iter1] = core_RMAMR(D,lambda,I,J);
end

M_hat = L + S + E;
error = norm(M_hat(:)-M(:))/norm(M(:));
disp(['Error: ' num2str(error)]);

%{
show_2dvideo(M,m,n);
show_2dvideo(M.*CM,m,n);
show_2dvideo(L,m,n);
show_2dvideo(S,m,n);
show_2dvideo(Z,m,n);
%}