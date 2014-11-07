%% [struct] = load_input(inputPath)
% video - struct
% 
function [video] = load_video_file(file)

% displog('Checking input file extension...');
input_extension = get_file_extension(file);

if(strcmp(input_extension,'mat'))
  % video2mat('dataset/demo.avi', 'dataset/demo.mat');
  % file = 'dataset/demo.mat';
  % displog('Loading data...');
  load(file);
else
  % displog('Reading movie file...');
  addpath('libs/mmread');
  video = mmread(file,[],[],false,true);
  [~,n] = size(video.frames);
  if(video.width > 320 || video.height > 320)
    width = round(video.width/2);
    height = round(video.height/2);
    for i = 1:n
      video.frames(i).cdata = imresize(rgb2gray(video.frames(i).cdata), [height width]);
    end
    video.width = width;
    video.height = height;
  else
    for i = 1:n
      video.frames(i).cdata = rgb2gray(video.frames(i).cdata);
    end
  end
  bytesize(video);
  % movie(video.frames);
  rmpath('libs/mmread');
  
  % For debug
  % show_video(video);
  
  %warning('on','all');
end
%displog('OK');

end