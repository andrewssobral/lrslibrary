function gui_show_results(files, axes, logMessageHandle)
  global stopFlag;
  global runningFlag;

  I_file = files.I_file;
  O_file = files.O_file;
  L_file = files.L_file;
  S_file = files.S_file;

  I_axes = axes.I_axes;
  O_axes = axes.O_axes;
  L_axes = axes.L_axes;
  S_axes = axes.S_axes;

  disp('Loading videos...');
  if(exist(I_file, 'file')) I_mov = load_video_file(I_file); I_mov = I_mov.frames; else I_mov = []; end
  if(exist(O_file, 'file')) O_mov = load_video_file(O_file); O_mov = O_mov.frames; else O_mov = []; end
  if(exist(L_file, 'file')) L_mov = load_video_file(L_file); L_mov = L_mov.frames; else L_mov = []; end
  if(exist(S_file, 'file')) S_mov = load_video_file(S_file); S_mov = S_mov.frames; else S_mov = []; end

  if(isempty(I_mov))
    warningMessage = sprintf('Warning: input file does not exist:\n%s', I_file);
    uiwait(msgbox(warningMessage));
    return;
  end

  [~,q] = size(I_mov);

  disp('Showing videos content...');
  runningFlag = true;

  for i = 1:q
    I = I_mov(i).cdata;
    imshow(I,'parent',I_axes);
    
    if(~isempty(O_mov))
      O = O_mov(i).cdata;
      imshow(O,'parent',O_axes);
    end
    
    if(~isempty(L_mov))
      L = L_mov(i).cdata;
      imshow(L,'parent',L_axes);
    end
    
    if(~isempty(S_mov))
      S = S_mov(i).cdata;
      imshow(S,'parent',S_axes);
    end
    
    set(logMessageHandle,'String',['Frame ' int2str(i)]);
    pause(0.01);

    if(stopFlag == true);
      break;
    end
  end
  runningFlag = false;

  disp('Finished!');
  stopFlag = false;
end
