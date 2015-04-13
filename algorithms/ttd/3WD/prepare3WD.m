function [D, A, E,O, numIter] = prepare3WD(confMap,detThresh,fileNames, numImages, params,destDir)
D = [] ;
for fileIndex = 1 : numImages
    currentImage = double(imread(fileNames{fileIndex}))/255;
    y   = currentImage(:);
    y = y / norm(y) ;
    D = [D y] ;
end

% Prepare object confidence matrix
M =  0.*D; 
for fileIndex = 1 : numImages-1
    ob = confMap(:,:,fileIndex);
    ob = ob(:);
    ind = find (ob>detThresh);
    mask = zeros(params.canonicalImageSize);
    mask = mask(:);
    mask(ind) = ob(ind);
    if ~isempty(ind)
        maskmax = max(mask);
        maskmin = min(mask);
        mask = (mask - maskmin)./(maskmax-maskmin);
        mm = reshape(mask,params.canonicalImageSize);
        mm2 = imdilate(mm,strel('ball',15,15));
        mm2max = max(mm2(:));
        mm2min = min(mm2(:));
        mm2n = (mm2 - mm2min)./(mm2max-mm2min);
         ind = find (mm2n>0);
         mm2n(ind) = mm2n(ind)+.5;
         subind = find (mm2n>1);
         mm2n(subind) = 1;
        mask = mm2n;
    end    
    M(:,fileIndex) = mask(:);
end

tau = params.tauc/sqrt(size(D,1));
lambda = params.lambdac/sqrt(size(D,1)) ; 
[A, E,O, numIter] = ThreeWayDec(D,M, tau,lambda, params.inner_tol, params.maxIter);
save(fullfile(destDir, 'final.mat'),'A','E','O','M') ;



