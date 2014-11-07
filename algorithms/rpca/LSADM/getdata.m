function opts = getdata(dataformat)

for i = 1:200
  imshow(reshape(images(:,i),144,176),[]); 
  pause(0.1);
end

if strcmp(dataformat,'surveillance-video-Hall')
    load Hall_airport_1000_1496_497_144_176_gray.mat;
    D = images(:,1:200); [n1,n2] = size(images); % imn1 = 144; imn2 = 176;
    opts.D = D; opts.mu = norm(D)/1.25;
    Xs = D; Ys = D;
    
elseif strcmp(dataformat,'surveillance-video-Campus-color')
%     load Campus_trees_1000_1993_994_128_160_gray.mat;
    load Campus_trees_1000_1996_997_128_160_R.mat; DR = images(:,181:500); 
    load Campus_trees_1000_1995_996_128_160_G.mat; DG = images(:,181:500);
    load Campus_trees_1000_1994_995_128_160_B.mat; DB = images(:,181:500);
    D = [DR, DG, DB];
    [n1,n2] = size(D); % imn1 = 128; imn2 = 160;
    opts.D = D; opts.mu = norm(D)/1.25;
    Xs = D; Ys = D;
end

opts.Xs = Xs;  opts.Ys = Ys;
opts.n1 = n1; opts.n2 = n2;
opts.sigma = 1e-6; opts.maxitr = 500; opts.rho = 1/sqrt(n1); 
opts.eta_mu = 2/3; opts.eta_sigma = 2/3;
opts.muf = 1e-6;
opts.sigmaf = 1e-6;
opts.epsilon = 1e-7;
opts.sv = 100;
