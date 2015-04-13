function gui_display_speed(n, handle)
  lrs_load_conf;
  
  fig_dir = fullfile(lrs_conf.lrs_dir,'figs');
  fig_name = fullfile(fig_dir,'speed0.png');
  
  switch n
    case 1
      fig_name = fullfile(fig_dir,'speed1.png');
    case 2
      fig_name = fullfile(fig_dir,'speed2.png');
    case 3
      fig_name = fullfile(fig_dir,'speed3.png');
    case 4
      fig_name = fullfile(fig_dir,'speed4.png');
    case 5
      fig_name = fullfile(fig_dir,'speed5.png');
  end
  
  img_speed = imread(fig_name,'BackgroundColor',[0.94 0.94 0.94]);
  imshow(img_speed,'parent',handle);
end