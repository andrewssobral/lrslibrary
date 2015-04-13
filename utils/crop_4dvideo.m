%%% [video_aux] = crop_4dvideo(video,r,s)
%
function [video_aux] = crop_4dvideo(video,r,s)
  if(s == -1)
    s = size(video,4);
  end
  
  k = 1;
  for i = r : s
    video_aux(:,:,:,k) = video(:,:,:,i);
    k = k + 1;
  end
end
