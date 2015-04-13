function fileNames = list_image_files(parameterString)

if nargin==0,
    %clc; clear;
    parameterString = '.';
end

temp = dir(parameterString);
fileNames = {};
[fileNames{1:length(temp),1}] = deal(temp.name);

% Find hidden files.
hiddenFiles = nan(size(fileNames));
for fileIndex = 1:length(fileNames),
    tempName = fileNames{fileIndex};
    hiddenFiles(fileIndex) = tempName(1) == '.';
end

% Find which entries in the listing are images.
imageFlags = nan(size(fileNames));
for fileIndex = 1:length(fileNames),
    [boo1, boo2, tempExtension] = fileparts(fileNames{fileIndex});
    tempExtension = lower(tempExtension);
    imageFlags(fileIndex) = ismember(tempExtension, {'.bmp', '.png', '.tif', '.jpg'});
end

% Find which entries in the listing are directories.
directoryFlags = {};
[directoryFlags{1:length(temp),1}] = deal(temp.isdir);
directoryFlags = cell2mat(directoryFlags);

tempFlags = (~hiddenFiles) & (~directoryFlags) & (imageFlags);
fileNames = fileNames(tempFlags);
