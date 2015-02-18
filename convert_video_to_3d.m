%%% [3dmatrix] = convert_video_to_3d(struct)
%
function [V] = convert_video_to_3d(video)
  for i = 1:video.nrFramesTotal
    frame = video.frames(i).cdata;
    if(size(frame,3) == 3)
      frame = rgb2gray(frame);
    end
    V(:,:,i) = frame;
  end
end
