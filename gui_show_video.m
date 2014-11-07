function gui_show_video(videoHandle, outputHandle, logMessageHandle)
  fullFileName = get(videoHandle,'String');
  gui_show_video_file(fullFileName, outputHandle, logMessageHandle);
end
