%% Download and install LRSLibrary
% https://github.com/andrewssobral/lrslibrary
run('../lrslibrary/lrs_setup');

%% Update the path
addpath('STOC-RPCA');

%% Run demo
folder = 'hyperspectral/';
folder_out = strcat(folder,'FG/');
folder_frames_rgb = strcat(folder,'imgRGB/');
folder_frames_tif = strcat(folder,'imgMS/');
d_frames_rgb = dir(strcat(folder_frames_rgb,'*.png'));
d_frames_tif = dir(strcat(folder_frames_tif,'*.tif'));
nframes = size(d_frames_rgb,1);
factor = 0.25; % input_size = [120 120];
k = 0;
for i = 1:20%nframes
  k = k + 1;
  disp(['frame: ' num2str(i)]);
  filename_rgb = strcat(folder_frames_rgb,d_frames_rgb(i).name);
  filename_tif = strcat(folder_frames_tif,d_frames_tif(i).name);
  
  input_rgb = im2double(imread(filename_rgb));
  input_tif = [];
  for j = 1:7, input_tif(:,:,j) = im2double(imread(filename_tif, j)); end
  
  input_rgb = imresize(input_rgb,factor);
  input_tif = imresize(input_tif,factor);
  %input_rgb = imresize(input_rgb,input_size);
  %input_tif = imresize(input_tif,input_size);
  
  % Only gray scale
    %frame = rgb2gray(input_rgb);
  % Only RGB values
    frame = input_rgb;
  % Only Multispectral bands
    %frame = input_tif;
  % Full: RGB + Multispectral bands
    %frame = cat(3,input_rgb,input_tif);
  
  % To see the frontal slices of the tensor
    %tensorlab.slice3(frame), colormap('gray'); 
  
  T = tensor(frame);
  if(k == 1) Tm = []; end
  [Tlowrank,Tsparse,Tmask,Tm] = OSTD(T,k,Tm);
  
  clf;
  subplot(1,4,1), imshow(input_rgb,[]);
  subplot(1,4,2), imshow(Tlowrank,[]);
  subplot(1,4,3), imshow(Tsparse,[]);
  subplot(1,4,4), imshow(Tmask,[]);
  pause(0.01);
  
  %%% Save outputs
  [~,filename_out,~] = fileparts(filename_rgb);
  filepath_out = strcat(folder_out,'fg',filename_out(3:end),'.png');
  %imwrite(Tmask,filepath_out); % uncomment for saving foreground masks
end % end frames
