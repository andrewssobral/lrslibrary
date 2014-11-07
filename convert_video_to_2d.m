%%% [2dmatrix] = convert_video_to_2d(struct)
%
function [I] = convert_video_to_2d(video)
  I = uint8(zeros(video.width*video.height,video.nrFramesTotal));
  
  for i = 1:video.nrFramesTotal
    frame = video.frames(i).cdata;
    I(:,i) = reshape(frame,[],1);
  end
end
