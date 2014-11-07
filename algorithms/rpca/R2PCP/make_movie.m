outputVideo = VideoWriter(fullfile('background-video.avi'));
outputVideo.FrameRate=12;
open(outputVideo);

A=reshape(A,[z1 z2 z3]); ov=imnorm(A);
% B=reshape(B,[z1 z2 z3]); ov=imnorm(B);
% Z=reshape(Z,[z1 z2 z3]); ov=imnorm(Z);

for ii = 1:size(ov,3)
    img = ov(:,:,ii);
    writeVideo(outputVideo,img);
end

close(outputVideo);



outputVideo = VideoWriter(fullfile('foreground-video.avi'));
outputVideo.FrameRate=12;
open(outputVideo);

% A=reshape(A,[z1 z2 z3]); ov=imnorm(A);
B=reshape(B,[z1 z2 z3]); ov=imnorm(B);
% Z=reshape(Z,[z1 z2 z3]); ov=imnorm(Z);

for ii = 1:size(ov,3)
    img = ov(:,:,ii);
    writeVideo(outputVideo,img);
end

close(outputVideo);



outputVideo = VideoWriter(fullfile('raw-video.avi'));
outputVideo.FrameRate=12;
open(outputVideo);

% A=reshape(A,[z1 z2 z3]); ov=imnorm(A);
% B=reshape(B,[z1 z2 z3]); ov=imnorm(B);
Z=reshape(Z,[z1 z2 z3]); ov=imnorm(Z);

for ii = 1:size(ov,3)
    img = ov(:,:,ii);
    writeVideo(outputVideo,img);
end

close(outputVideo);


disp('videos written.')
clear outputVideo img ov
