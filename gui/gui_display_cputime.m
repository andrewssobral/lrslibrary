function gui_display_cputime(n, handle)
  lrs_load_conf;
  
  fig_dir = fullfile(lrs_conf.lrs_dir,'figs');
  fig_name = fullfile(fig_dir,'time0.png');
  
  switch n
    case 1
      fig_name = fullfile(fig_dir,'time1.png');
    case 2
      fig_name = fullfile(fig_dir,'time2.png');
    case 3
      fig_name = fullfile(fig_dir,'time3.png');
    case 4
      fig_name = fullfile(fig_dir,'time4.png');
    case 5
      fig_name = fullfile(fig_dir,'time5.png');
  end
  
  img_time = imread(fig_name,'BackgroundColor',[0.94 0.94 0.94]);
  imshow(img_time,'parent',handle);
end