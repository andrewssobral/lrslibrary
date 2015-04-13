function saveResults(destDir, canonicalImageSize)
% alignment results
load(fullfile(destDir, 'final.mat'),'A','E','O','M') ;

%% display
mkdir([destDir '\Err']);
mkdir([destDir '\Res']);
mkdir([destDir '\O']);
FACT = 255;
for i=1:size(A,2)
    im = reshape(A(:,i), canonicalImageSize);
    imMin = min(im(:));
    imMax = max(im(:));
    im = (im - imMin)./(imMax - imMin);
    im = uint8(255*im);
    imwrite(im,[destDir '\Res\' sprintf('Image_%.4d.bmp',i)],'bmp');
    
    % E
    im = reshape(E(:,i), canonicalImageSize);
    im = abs(im) ./ max(max(abs(im)));
    im = uint8(FACT.*abs(im));
    imwrite(im,[destDir '\Err\' sprintf('Er_%.4d.bmp',i)],'bmp');
    
    % O
    im = reshape(O(:,i), canonicalImageSize);
    im = abs(im) ./ max(max(abs(im)));
    im = uint8(FACT.*abs(im));
    imwrite(im,[destDir '\O\' sprintf('O_%.4d.bmp',i)],'bmp');
end;



