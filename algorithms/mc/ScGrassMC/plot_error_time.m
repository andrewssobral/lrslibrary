% meth_objs: cell of column vectors (of different lengths)
% meth_names: method names
function plot_error_time(fname, meth_objs, meth_times, meth_names, titles, scale, disptime)

if ~exist('disptime', 'var')
  disptime = 0;
  for meth_id = 1:length(meth_objs)
    disptime = max(disptime,length(meth_objs{meth_id}));
  end
end

num_meth = length(meth_objs);
line_styles = {'ks-.', 'rx-', 'bv--', 'go-', 'y*-', 'c-.'};
line_styles = {'s-', 'x-', 'v-', 'o-', '*-', 'v-','+-','d-','^-','<-','>-','p-','h-'};
colors = distinguishable_colors(num_meth,[1, 1, 1]);

maxtime = 0;
for i=1:num_meth
  maxtime = max(maxtime, max(meth_times{i}));
end

for j=1:size(meth_objs{1},2)
  h = figure;
  hold on;
  for i = 1:num_meth
    times = meth_times{i}(1:min(size(meth_times{i},j),size(meth_times{i},j)));
    data  = meth_objs{i}(1:min(size(meth_times{i},j),size(meth_times{i},j)));
    if strcmp(scale,'log')
      data = log(data);
      %plot(times, data(:,j),line_styles{i},'LineWidth',1,'Color',colors(i,:));
    end
    plot(times, data(:,j),line_styles{i},'LineWidth',1,'Color',colors(i,:));
  end
  legend(meth_names,'Location','Best');
  hold off;
end

xlim([0,min(disptime,maxtime)]);
xlabel('Time [s]');
ylabel('RMSE');
if strcmp(scale,'log')
  ylabel('RMSE (log-scale)');
end
title(titles);

if ~isempty(fname)
  print(h,'-dpsc',fname);
end

