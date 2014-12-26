function U = initializeSubspace(inputFolder, imgDir, param)
    % initial lize the  subspace by svd    
   
    trainSize = min(length(imgDir), param.trainSize);
    sampleSize = min(param.sampleSize, trainSize);

    randIndex =   randperm(trainSize);
    
    imageIndex =  randIndex(1:sampleSize);
    
    %imageIndex =  randperm(trainSize, sampleSize);
    
    U0=[];
    for i=1:sampleSize
        index = imageIndex(i);
        imgPath = [inputFolder, imgDir(index).name];
        im = imread(imgPath);
        im = im(:);
        U0 = [U0 im];
    end
    
    [U, D] = svds(double(U0), param.rank);
end