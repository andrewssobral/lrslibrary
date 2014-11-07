function [tr,vd,ts]=RandPartDB(gnd,trNum)
% partition dataset in three subsets
% gnd, class labels
% trNum, feature number selected individually
% copyright @ Guan Naiyang
% modified time:
% 18/4/2010

lb = unique(gnd);
n = length(gnd);
nn = length(lb);
tr = [];
ts = [];
for j = 1:nn
    p_index = find(gnd==lb(j));
    index = randperm(length(p_index));
    tr = [tr, p_index(index(1:trNum))'];
    ts = [ts, p_index(index(trNum+1:(2*trNum)))'];
end
vd = 1:n;
vd([tr,ts]) = [];

return;