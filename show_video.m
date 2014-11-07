%%% void show_video(struct)
%
function show_video(video)
  for i = 1 : video.nrFramesTotal
    frame = video.frames(i).cdata;
    imshow(frame,[],'InitialMagnification','fit');
    disp(i);
    pause(0.01);
  end
end
