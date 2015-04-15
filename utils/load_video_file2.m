%% [3dmatrix/struct] = load_video_file2(string,struct)
%
function [V,obj] = load_video_file2(filename,params)
  color = 'rgb';
  nFrames = -1;
  format = 'default';
  debug = 0;
  shape = 'S';
  
  if(isfield(params,'color')) color = params.color; end
  if(isfield(params,'nFrames')) nFrames = params.nFrames; end
  if(isfield(params,'format')) format = params.format; end
  if(isfield(params,'debug')) debug = params.debug; end
  if(isfield(params,'shape')) shape = params.shape; end
  
  obj = VideoReader(filename);
  
  if(nFrames <= 0) nFrames = obj.NumberOfFrames; end
  if(nFrames > obj.NumberOfFrames) nFrames = obj.NumberOfFrames; end
  
  k = 1;
  for i = 1:nFrames
    frame = read(obj, i);
    
    if(isfield(params,'rs')) frame = imresize(frame,params.rs); end
    
    if(strcmp(color,'gray')) frame = rgb2gray(frame); end
    if(strcmp(color,'hsv')) frame = rgb2hsv(frame); end
    
    if(strcmp(format,'double')) frame = im2double(frame); end
    
    if(strcmp(shape,'S')) V{k} = frame; end
    if(strcmp(shape,'2D')) V(:,k) = reshape(frame,[],1); end % only for gray
    if(strcmp(shape,'3D')) V(:,:,k) = frame; end % only for gray
    if(strcmp(shape,'4D')) V(:,:,:,k) = frame; end
    
    if(debug == 1)
      imshow(frame,[],'InitialMagnification','fit');
      disp(i);
      pause(0.1);
    end
    
    k = k + 1;
  end
end
