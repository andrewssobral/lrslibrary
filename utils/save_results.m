%%% void save_results(struct,string)
%
function save_results(movobj,outFile)
  if(strcmp(get_file_extension(outFile),'mat'))
    disp('Saving in mat file');
    save(outFile,'movobj');
  else  
    disp('Saving results in movie file');

    if(~isempty(movobj.L))
      L_file = gen_file_name(outFile,'_L');
      disp(['Saving low rank results at: ' L_file]);
      v = VideoWriter(L_file,'Uncompressed AVI');
      open(v);
      writeVideo(v,movobj.L);
      close(v);
      %movie2avi(movobj.L, L_file, 'compression', 'None');
      clear L_file;
    end

    if(~isempty(movobj.S))
      S_file = gen_file_name(outFile,'_S');
      disp(['Saving sparse results at: ' S_file]);
      v = VideoWriter(S_file,'Uncompressed AVI');
      open(v);
      writeVideo(v,movobj.S);
      close(v);
      %movie2avi(movobj.S, S_file, 'compression', 'None');
      clear S_file;
    end

    if(~isempty(movobj.O))
      disp(['Saving foreground result at: ' outFile]);
      v = VideoWriter(outFile,'Uncompressed AVI');
      open(v);
      writeVideo(v,movobj.O);
      close(v);
      %movie2avi(movobj.O, outFile, 'compression', 'None');
    end
  end
end
