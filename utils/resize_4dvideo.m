%%% [video_aux] = resize_4dvideo(video,double)
% f : [0,...,1,...,N]
function [video_aux] = resize_4dvideo(video,f)
  [~, ~, ~, nFrames] = size(video);
  for i = 1 : nFrames
    frame = video(:,:,:,i);
    frame = imresize(frame,f);
    video_aux(:,:,:,i) = frame;
  end
end
