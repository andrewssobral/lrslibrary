function newLabel=GetLabel(trnData,tstData,trnLabel)

% Shared by algorithms.
% copyright @ guan naiyang

dist=EuDist2(tstData',trnData',0);
[junk, sortedIdx]=sort(dist,2);
newLabel=trnLabel(sortedIdx);
clear dist;

return;