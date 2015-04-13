function [fileNames, numImages] = get_training_images(imagePath)
fileNames = {};
imageIndex = 0;
userDirectoryContents = list_image_files(imagePath);
if isempty(userDirectoryContents)
    error(['No image files were found! Check your paths; there should be images in ' fullfile(rootPath, trainingDatabaseName)]);
end
for fileIndex = 1:length(userDirectoryContents),
    imageName = userDirectoryContents{fileIndex};
    disp(['Using image file ' imageName '...']);
    
    imageIndex = imageIndex+1;
    
    imageFileName = fullfile(imagePath, imageName);
    fileNames{imageIndex} = imageFileName;
end
numImages = length(userDirectoryContents);
