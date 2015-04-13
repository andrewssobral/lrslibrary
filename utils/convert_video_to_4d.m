%%% [4dmatrix] = convert_video_to_4d(struct)
%
function [V] = convert_video_to_4d(video)
  for i = 1:video.nrFramesTotal
    frame = video.frames(i).cdata;
    V(:,:,:,i) = frame;
  end
end
