clear, clc;
inFile='dataset/demo.avi';
video = load_video_file(inFile);
opts.rows = video.height;
opts.cols = video.width;
M = im2double(convert_video_to_2d(video));
% imagesc(M); colormap('gray');

%%
regularizer = 0.1;
param.superpixel.slicParam = [10, regularizer; 20, regularizer; 40,regularizer; 80, regularizer;];
%  param.superpixel.slicParam = [5,regularizer;  10, regularizer; 20, regularizer;  40,regularizer ; ];
param.maxIter = 200;
param.tol = 1e-4;
param.lambda = 0.33;    %   bestAccu = 0.9762, lambda = 0.33; threshold = 8e-6;
param.rank = 2; % subspace dimension
param.eta = 1e-3;   % stepsize for subspace updating
%param.sampleSize = 20; % number of sample frames to learn the background
%param.trainSize = 24; % sample from the fist 240 frames to learn the background
%param.startIndex = 1;   %  which frame to start
%param.randomStart = true;  
param.admm.I = speye(size(M,1)*3);

%%
%U = orth(randn(size(M,1)*3, param.rank));
[U, ~] = svds([M(:,1:25);M(:,1:25);M(:,1:25)], param.rank);

%%
run('C:/GitHub/lrslibrary/libs/vlfeat-0.9.19/toolbox/vl_setup')

%%
for i = 1:size(M,2)
  disp(i);
  I = M(:,i);
  im = im2uint8(reshape(I,[opts.rows opts.cols]));
  im = repmat(im,[1 1 3]);
  
  param.admm.G = getGroupSuperColor27(im, param);
  param.admm.Z = sparse(size(M,1)*3, size(param.admm.G, 2), 0);
  param.admm.Y = param.admm.Z;
  
  v = double(im(:))/255.0;
  
  [x, w] = solveWXADMM(U, v, param);

  pL = U*w; % low-rank
  pS = pL-v; % sparse
  %E = v;

  imgL = reshape(pL, size(im));
  imgS = reshape(pS, size(im));

  vL = imgL(:,:,1);
  vS = imgS(:,:,1);

  vL = vL(:);
  vS = vS(:);

  L(:,i) = vL;
  S(:,i) = vS;
  
  %L = U*w; % low-rank
  %E = v;
  
  %fg = L-v;
  %fg(abs(x)<=8e-10) = 0;
  %fg(abs(x)>8e-10) = 1;
  
  %subplot(1,3,1), imshow(im,[]);
  %subplot(1,3,2), imshow(reshape(L, size(im)),[]);
  %subplot(1,3,3), imshow(reshape(fg, size(im)),[]);
  
  residual = param.lambda*(U*w + x - v);
  U = updateSubspace(U, residual, w, param);
  
  %pause(.01);
end
