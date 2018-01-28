%% [struct] = load_input(inputPath)
% video - struct
% 
function [video] = load_video_file(file)

% displog('Checking input file extension...');
input_extension = get_file_extension(file);
%constraint = 0;

if(strcmp(input_extension,'mat'))
  % video2mat('dataset/demo.avi', 'dataset/demo.mat');
  % file = 'dataset/demo.mat';
  % displog('Loading data...');
  load(file);
else
  % displog('Reading movie file...');
  
  %video = mmread(file,[],[],false,true);
  xyloObj = VideoReader(file);
  vidWidth = xyloObj.Width;
  vidHeight = xyloObj.Height;
  k = 0;
	while hasFrame(xyloObj)
    k = k+1;
    video.frames(k).cdata = readFrame(xyloObj);
  end
  video.nrFramesTotal = k;
  video.width  = vidWidth;
  video.height = vidHeight;

%   if(constraint)
%     if(video.width > 320 || video.height > 320)
%       width = round(video.width/2);
%       height = round(video.height/2);
%       [~,n] = size(video.frames);
%       for i = 1:n
%         video.frames(i).cdata = imresize(video.frames(i).cdata, [height width]);
%       end
%       video.width = width;
%       video.height = height;
%     end
%   end
  
  bytesize(video);
  % movie(video.frames);
  
  % For debug
  % show_video(video);
  
  %warning('on','all');
end
%displog('OK');

end