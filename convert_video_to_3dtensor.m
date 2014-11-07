%%% [tensor] = convert_video_to_3dtensor(struct)
%
function [T] = convert_video_to_3dtensor(video)
  V = convert_video_to_3d(video);
  T = tensor(V);
end
