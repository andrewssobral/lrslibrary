function gui_show_video_file(fullFileName, outputHandle, logMessageHandle)
  global stopFlag;
  global runningFlag;

  disp(['Input video: ' fullFileName]);

  if exist(fullFileName, 'file')
    % File exists.  Do stuff....
    disp('Loading video...');
    video = load_video_file(fullFileName);
    mov = video.frames;
    [p,q] = size(mov);

    disp('Showing video content...');
    runningFlag = true;
    for i = 1:q
      B = mov(i).cdata;
      imshow(B,'parent',outputHandle);
      set(logMessageHandle,'String',['Frame ' int2str(i)]);
      pause(0.01);

      if(stopFlag == true);
          break;
      end
    end
    runningFlag = false;

    disp('Finished!');
    stopFlag = false;
  else
    % File does not exist.
    warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
    uiwait(msgbox(warningMessage));
  end
end

