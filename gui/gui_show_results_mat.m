function gui_show_results_mat(I_file, movs, axes, logMessageHandle)
  global stopFlag;
  global runningFlag;

  disp('Loading videos...');

  if(exist(I_file, 'file')) I_mov = load_video_file(I_file); I_mov = I_mov.frames; else I_mov = []; end

  if(isempty(I_mov))
    warningMessage = sprintf('Warning: input file does not exist:\n%s', I_file);
    uiwait(msgbox(warningMessage));
    return;
  end

  O_mov = movs.O_mov;
  L_mov = movs.L_mov;
  S_mov = movs.S_mov;

  I_axes = axes.I_axes;
  O_axes = axes.O_axes;
  L_axes = axes.L_axes;
  S_axes = axes.S_axes;

  [~, q] = size(I_mov);

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
