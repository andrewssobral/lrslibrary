function gosus(param)

     
    dataPath =  param.input.dataPath;

    dataset = param.input.dataset;

    inputFolder = [dataPath, filesep, dataset, filesep];

    dirString = param.input.dirString;

    resultPath = param.output.resultPath;
 
    outputFolder = [resultPath, filesep, dataset, filesep];
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
 
    imgDir = dir([inputFolder,filesep, dirString]);
    threshold = param.output.threshold;
    
    
     
    imgPath = [inputFolder, imgDir(1).name];
    im = imread(imgPath);
    [rows,cols, channel] = size(im);
    param.admm.I = speye(rows*cols*channel);    

    if  param.randomStart 
        U =  orth(randn(rows*cols*channel, param.rank));
    else
        U = initializeSubspace(inputFolder, imgDir, param);
    end

    startIndex = param.startIndex;        
    
    if param.showResult
            hOrig = subplot(1,3,1);
            set(gca,'nextplot','replacechildren'); 
            title('Original Image');
            hBg = subplot(1,3,2);
            set(gca,'nextplot','replacechildren'); 
            title('Background');            
            hFg = subplot(1,3,3);
            set(gca,'nextplot','replacechildren'); 
            title('Foreground');
    end
    
 
    for i=startIndex:  length(imgDir)
     
        objName = imgDir(i).name
        
        imgPath = [inputFolder, imgDir(i).name];
        im = imread(imgPath);
        
        param.admm.G = getGroupSuperColor27(im, param);  
        param.admm.Z =   sparse( rows*cols*channel,  size(param.admm.G, 2), 0 );
        param.admm.Y = param.admm.Z;     

        v =double(im(:))/255.0;
  
        [  x, w ] = solveWXADMM( U, v, param );
 
        if param.saveResult  % save segmented foreground and background
            saveResult;
        end        
        
        if param.showResult
            showResult;
        end

%         residual = U*w + x - v;
        residual = param.lambda*( U*w + x - v);
        
        U = updateSubspace(U, residual, w, param);     

    end
    
    
    function showResult()
        
        imOrig = reshape(v, size(im)); 
        axes(hOrig);
        imshow( imresize(imOrig,1) );
       

        bg =U*w;       
        imBG = reshape(bg, size(im));
        axes(hBg);
        imshow( imresize(imBG,1) );
       
         
        fg = v;
        fg(abs(x)<=threshold) = 0;
        imFG = reshape(fg, size(im)); 
        axes(hFg);
        imshow( imresize(imFG,1) );
        
   
    end

   function saveResult()
        
       
        imOrig = reshape(v, size(im)); 
        imwrite( imresize(imOrig,1) , [outputFolder, objName(1:end-4), '_orig.jpg'], 'jpg');

        fg = v;
        fg(abs(x)<=threshold) = 0;
        imFG = reshape(fg, size(im)); 
        imwrite( imresize(imFG,1) , [outputFolder, objName(1:end-4), '_fg.jpg'], 'jpg');

        bg =U*w;
        imBG = reshape(bg, size(im));
        imwrite( imresize(imBG,1) , [outputFolder, objName(1:end-4), '_bg.jpg'], 'jpg');
      
    end
end