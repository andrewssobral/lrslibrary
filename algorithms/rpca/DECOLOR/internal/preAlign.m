function [ImTrans,tau] = preAlign(ImData)
% for pre-alignment
numFrame = size(ImData,3);
ImTrans = ImData;
IDcenter = round(numFrame/2);
tau = zeros(6,numFrame);
numLevel = ceil(log2(max([size(ImData,1),size(ImData,2)])/50))+1;
for i = IDcenter-1:-1:1
    disp(['frame: ',num2str(i)]);
    [ImTrans(:,:,i),tau(:,i)] = regMGNC(ImTrans(:,:,i+1),ImData(:,:,i),tau(:,i+1),numLevel);
end
for i = IDcenter+1:+1:numFrame
    disp(['frame: ',num2str(i)]);
    [ImTrans(:,:,i),tau(:,i)] = regMGNC(ImTrans(:,:,i-1),ImData(:,:,i),tau(:,i-1),numLevel);
end
end
