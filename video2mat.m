%%% void = video2mat(string, string)
% input - string
% output - string
%
% i.e: video2mat('dataset/video.avi', 'dataset/video.mat')
%
function video2mat(input, output)
  video = load_video_file(input);
  save(output,'video');
  disp('Saved!');
end
