%% void displog(string)
% msg - string
%
function displog(msg)
  disp([datestr(now, 'HH:MM:SS') ' ' msg]);
end
