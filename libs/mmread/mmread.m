function [video, audio] = mmread(filename, frames, time, disableVideo, disableAudio, matlabCommand, trySeeking, useFFGRAB)
% [video, audio] = mmread(filename, frames, time, disableVideo, 
%                       disableAudio, matlabCommand, trySeeking, useFFGRAB)
% mmread reads virtually any media file.  It now uses AVbin and FFmpeg to 
% capture the data, this includes URLs.  The code supports all major OSs
% and architectures that Matlab runs on.
%
% INPUT
% filename      input file to read (mpg, avi, wmv, asf, wav, mp3, gif, ...)
% frames        specifies which video frames to capture, default [] for all
%               or to specify time
% time          [startTime stopTime], default [] for all
% disableVideo  disables ALL video capturing, to save memory or time
% disableAudio  disables ALL audio capturing, to save memory or time
% matlabCommand Do not return the video structure, but call the function
%               specified by matlabCommand.  The function definition must
%               match that of processFrame.m.  See processFrame.m for more
%               information.
% trySeeking    [true] setting this to false makes the code slower but more
%               precise.  If the first several frames are distorted or
%               timing information isn't accurate, set this to false.
% useFFGRAB     [true] Use the new version of mmread, which uses ffmpeg.
%               However, if an audio or video stream can't be read AND you 
%               are running Windows try setting this to false (old version).
%
% OUTPUT
% video is a struct with the following fields:
%   width           width of the video frames
%   height          height of the video frames
%   rate            the frame rate of the video, if it can't be determined
%                   it will be 1.
%   nrFramesTotal   the total number of frames in the movie regardless of
%                   how many were captured.  Unfortunately, this can not
%                   always be determined.  If it is negative then it
%                   is an estimate based upon the duration and rate
%                   (normally accurate to within .1%).   It can be 0,
%                   in which case it could not be determined at all.  If it
%                   is a possitive number then it should always be accurate.
%   totalDuration   the total length of the video in seconds.
%   frames          a struct array with the following fields:
%       cdata       [height X width X 3] uint8 matricies
%       colormap    always empty
%   times           the corresponding time stamps for the frames (in msec)
%   skippedFrames   some codecs (not mmread) will skip duplicate frames
%                   (i.e. identical to the previous) in fixed frame rate
%                   movies to save space and time.  These skipped frames
%                   can be detected by looking for jumps in the "times"
%                   field.  This field will be true when frames are
%                   skipped.
%
% audio is a struct with the following fields:
%   nrChannels      the number of channels in the audio stream (1 or 2)
%   rate            sampling rate of the audio, ex. 44100.  If it can't be
%                   determined then it will be 1.
%   bits            bit depth of the samples (8 or 16)
%   data            the real data of the whole audio stream.  This can be
%                   played using wavplay.  If time ranges are specified,
%                   the length of the data may not correspond to the total
%                   time.  This normally happens with movies.  The issue is
%                   that the start of the audio stream is generally counted
%                   at the END of the first frame.  So, time is shifted by
%                   1/framerate.
%   nrFramesTotal   Audio comes in packets or frames when captured, the
%                   division of the audio into frames may or may not make
%                   sense.
%   totalDuration   the total length of the audio in seconds.
%   frames          cell array of uint8s.  Probably not of great use.
%   times           the corresponding time stamps for the frames (in milliseconds)
%
% If there is no video or audio stream the corresponding structure will be
% empty.
%
% Specifying frames does not effect audio capturing.  If you want only a
% subsection of the audio use the 3rd parameter "time".  Specifying time
% effects both audio and video.  Time is specified in seconds (subsecond
% resolution is supported with fractional numbers ex. 1.125), starting at 0.
% Time is defined as startTime (inclusive) to stopTime (exclusive), or
% using set notation [startTime stopTime).
%
% If there are multiple video or audio streams, then the structure will be
% of length > 1.  For example: audio(1).data and audio(2).data.
%
% EXAMPLES
% [video, audio] = mmread('chimes.wav'); % read whole wav file
% wavplay(audio.data,audio.rate);
%
% video = mmread('mymovie.mpg'); % read whole movie
% movie(video.frames);
%
% video = mmread('mymovie.mpg',1:10); %get only the first 10 frames
%
% video = mmread('mymovie.mpg',[],[0 3.5]); %read the first 3.5 seconds of the video
%
% [video, audio] = mmread('chimes.wav',[],[0 0.25]); %read the first 0.25 seconds of the wav
% [video, audio] = mmread('chimes.wav',[],[0.25 0.5]); %read 0.25 to 0.5 seconds of the wav, there is no overlap with the previous example.
%
% read a movie directly from a URL  
% video = mmread('http://www.nature.com/neuro/journal/v9/n4/extref/nn1660-S8.avi');  
%
% video = mmread('mymovie.mpg',[],[],false,true); %read all frames, disable audio
%
% mmread('mymovie.mpg',[],[],false,false,'processFrame'); %Use inline processing for all frames in a movie using the function processFrame.m
%
% Copyright 2008 Micah Richert
% 
% This file is part of mmread.

if nargin < 8
    useFFGRAB = true;
    if nargin < 7
        trySeeking = true;
        if nargin < 6
            matlabCommand = '';
            if nargin < 5
                disableAudio = nargout < 2;
                if nargin < 4
                    disableVideo = false;
                    if nargin < 3
                        time = [];
                        if nargin < 2
                            frames = [];
                        end
                    end
                end
            end
        end
    end
end

if useFFGRAB
    currentdir = pwd;
    try
        if ~ispc
            cd(fileparts(mfilename('fullpath'))); % FFGrab searches for AVbin in the current directory
        end

	fmt = '';
	if iscell(filename)
            if length(filename) ~= 2
                error('If you are specifying filename and format, they must be in a cell array of lenght 2.');
            end
            fmt = filename{2};
            filename = filename{1};
	end

        FFGrab('build',filename,fmt,double(disableVideo),double(disableAudio),double(trySeeking));
        
        if (isempty(time))
            FFGrab('setFrames',frames);
        else
            if (numel(time) ~= 2)
                error('time must be a vector of length 2: [startTime stopTime]');
            end
            FFGrab('setTime',time(1),time(2));
        end
        FFGrab('setMatlabCommand',matlabCommand);

        try
            FFGrab('doCapture');
        catch
            err = lasterror;
            if (~strcmp(err.identifier,'processFrame:STOP'))
                rethrow(err);
            end
        end

        [nrVideoStreams, nrAudioStreams] = FFGrab('getCaptureInfo');

        video = struct('width',{},'height',{},'nrFramesTotal',{},'frames',{});
        audio = struct('nrChannels',{},'rate',{},'bits',{},'nrFramesTotal',{},'data',{},'frames',{});

        warned = false;

        % we can only get the video frames if we don't process a matlabCommand
        if strcmp(matlabCommand,'')
            % loop through getting all of the video data from each stream
            for i=1:nrVideoStreams
                [width, height, rate, nrFramesCaptured, nrFramesTotal, totalDuration] = FFGrab('getVideoInfo',i-1);
                video(i).width = width;
                video(i).height = height;
                video(i).rate = rate;
                video(i).nrFramesTotal = nrFramesTotal;
                video(i).totalDuration = totalDuration;
                video(i).frames = struct('cdata',cell(1,nrFramesCaptured),'colormap',cell(1,nrFramesCaptured));
                video(i).times = zeros(size(video(i).frames));
                video(i).skippedFrames = [];

                if (nrFramesTotal > 0 && any(frames > nrFramesTotal))
                    warning('mmread:general',['Frame(s) ' num2str(frames(frames>nrFramesTotal)) ' exceed the number of frames in the movie.']);
                end

                scanline = ceil(width*3/4)*4; % the scanline size must be a multiple of 4.

                for f=1:nrFramesCaptured
                    [data, time] = FFGrab('getVideoFrame',i-1,f-1);

                    if any(size(data) == 0)
                        warning('mmread:getVideoFrame',['Frame ' num2str(f) ' could not be decoded']);
                    else
                        % the data ordering is wrong for matlab images, so permute it
                        data = permute(reshape(data, 3, width, height),[3 2 1]);
                        video(i).frames(f).cdata = data;
                        video(i).times(f) = time;
                    end
                end

                framerate = (max(video(i).times)-min(video(i).times))/nrFramesCaptured;
                if framerate > 0
                    video(i).skippedFrames = any(diff(video(i).times)>framerate*1.8) & abs(mean(diff(video(i).times))-framerate)/framerate<0.05;
                end

                % if frames are specified then make sure that the order is the same
                if (~isempty(frames) && nrFramesCaptured > 0)
                    [uniqueFrames, dummy, frameOrder] = unique(frames);
                    if (length(uniqueFrames) > nrFramesCaptured)
                        warning('mmread:general','Not all frames specified were captured.  Returning what was captured, but order may be different than specified.');
                        remainingFrames = frames(frames<=uniqueFrames(nrFramesCaptured));
                        [dummy, dummy, frameOrder] = unique(remainingFrames);
                    end

                    video(i).frames = video(i).frames(frameOrder);
                    video(i).times = video(i).times(frameOrder);
                end
            end
        end

        % loop through getting all of the audio data from each stream
        for i=1:nrAudioStreams
            [nrChannels, rate, bits, nrFramesCaptured, nrFramesTotal, subtype, totalDuration] = FFGrab('getAudioInfo',i-1);
            audio(i).nrChannels = nrChannels;
            audio(i).rate = rate;
            audio(i).bits = bits;
            audio(i).nrFramesTotal = nrFramesTotal;
            audio(i).totalDuration = totalDuration;
            audio(i).frames = cell(1,nrFramesCaptured);
            for f=1:nrFramesCaptured
                [data, time] = FFGrab('getAudioFrame',i-1,f-1);
                audio(i).frames{f} = data;
                audio(i).times(f) = time;
            end
            % combine the data across frames
            d = double(cat(1,audio(i).frames{:}));

            % rescale the data so that it is between -1.0 and 1.0
            if (subtype==0)
                %PCM formated data...
                switch (bits)
                    case {4, 8}
                        d = (d-2^(bits-1))/2^(bits-1);
                    case {16, 24, 32}
                        d = d/2^(bits-1);
                end
            elseif (subtype==1)
                if (bits == 32)
                    %IEEE FLOAT formated data...
                    if (max(d) > 1 | min(d) < -1)
                        % there are two float formats one that is already -1 to 1
                        % and the there is between -2^15 to 2^15
                        d = d / 2^15;
                    end
                else
                    warning('Audio data format not recognized/supported, it probably is going to be useless.');
                end
            else
                warning('Audio data format not recognized/supported, it probably is going to be useless.');
            end

            % reshape the data so that it is nrChannels x Samples.  This should be the same output as wavread.
            audio(i).data = reshape(d,nrChannels,length(d)/nrChannels)';
        end

        FFGrab('cleanUp');
    catch
        try
            err = lasterror;
            if ~isempty(strfind(err.message,'libavbin'))
                switch mexext
                    case 'mexa64'
%                        if strfind(err.message,'No such file')
                            if exist('libavbin.so.64','file')
                                if exist('libavbin.so','file')
                                    d64=dir('libavbin.so.64');
                                    d=dir('libavbin.so');
                                    if d.bytes ~= d64.bytes
                                        R=input('libavbin.so is installed but seems to be the 32bit version.\nShall I correct this (no admin required)? [Y/n]','s');
                                        if ~isequal(R,'n') & ~isequal(R,'N')
                                            if copyfile('libavbin.so.64','libavbin.so','f')
                                                error('libavbin.so installed.  You may need to restart Matlab for mmread to function.');
                                            else
                                                error('libavbin.so failed to install:  couldn''t write to libavbin.so');
                                            end
                                        end
                                    end
                                else
                                    R=input('libavbin.so needs to be installed.\nShall I install this for you (no admin required)? [Y/n]','s');
                                    if ~isequal(R,'n') && ~isequal(R,'N')
                                        if copyfile('libavbin.so.64','libavbin.so','f')
                                            error('libavbin.so installed.  You may need to restart Matlab for mmread to function.');
                                        else
                                            error('libavbin.so failed to install:  couldn''t write to libavbin.so');
                                        end
                                    end
                                end
                            end
%                        end
                    case 'mexglx'
%                        if strfind(err.message,'No such file')
                            if exist('libavbin.so.32','file')
                                if exist('libavbin.so','file')
                                    d32=dir('libavbin.so.32');
                                    d=dir('libavbin.so');
                                    if d.bytes ~= d32.bytes
                                        R=input('libavbin.so is installed but seems to be the 64bit version.\nShall I correct this (no admin required)? [Y/n]','s');
                                        if ~isequal(R,'n') && ~isequal(R,'N')
                                            if copyfile('libavbin.so.32','libavbin.so','f')
                                                error('libavbin.so installed.  You may need to restart Matlab for mmread to function.');
                                            else
                                                error('libavbin.so failed to install:  couldn''t write to libavbin.so');
                                            end
                                        end
                                    end
                                else
                                    R=input('libavbin.so needs to be installed.\nShall I install this for you (no admin required)? [Y/n]','s');
                                    if ~isequal(R,'n') && ~isequal(R,'N')
                                        if copyfile('libavbin.so.32','libavbin.so','f')
                                            error('libavbin.so installed.  You may need to restart Matlab for mmread to function.');
                                        else
                                            error('libavbin.so failed to install:  couldn''t write to libavbin.so');
                                        end
                                    end
                                end
                            end
%                        end
                    case {'mexmac', 'mexmaci', 'mexmaci64'}
                        R=input('libavbin.dylib needs to be installed.\nShall I install this for you (no admin required)? [Y/n]','s');
                        if ~isequal(R,'n') && ~isequal(R,'N')
                            if ~exist('~/lib','dir')
                                if ~mkdir('~/lib')
                                    error('A lib directory doesn''t exist in your home directory, and one couldn''t be created.');
                                end
                            end
                            if copyfile('libavbin.dylib','~/lib/libavbin.dylib','f')
                                error('libavbin.dylib installed.  You may need to restart Matlab for mmread to function.');
                            else
                                error('libavbin.dylib failed to install:  couldn''t write to libavbin.dylib');
                            end
                        end
                end
            end
            try
                FFGrab('cleanUp');
            catch
            end
            rethrow(err);
        catch
            if ~strmatch(computer,'PCWIN')
                cd(currentdir);
            end
            rethrow(lasterror);
        end
    end
    if ~ispc
        cd(currentdir);
    end
else
    try
        mexDDGrab('buildGraph',filename);
        if (isempty(time))
            mexDDGrab('setFrames',frames);
        else
            if (numel(time) ~= 2)
                error('time must be a vector of length 2: [startTime stopTime]');
            end
            mexDDGrab('setTime',time(1),time(2));
        end
        if (disableVideo)
            mexDDGrab('disableVideo');
        end;
        if (disableAudio | nargout < 2)
            mexDDGrab('disableAudio');
        end;
        mexDDGrab('setMatlabCommand',matlabCommand);

        mexDDGrab('setTrySeeking',double(trySeeking));

        try
            mexDDGrab('doCapture');
        catch
            err = lasterror;
            if (~strcmp(err.identifier,'processFrame:STOP'))
                rethrow(err);
            end
        end

        [nrVideoStreams, nrAudioStreams] = mexDDGrab('getCaptureInfo');

        video = struct('width',{},'height',{},'nrFramesTotal',{},'frames',{});
        audio = struct('nrChannels',{},'rate',{},'bits',{},'nrFramesTotal',{},'data',{},'frames',{});

        warned = false;

        % we can only get the video frames if we don't process a matlabCommand
        if strcmp(matlabCommand,'')
            % loop through getting all of the video data from each stream
            for i=1:nrVideoStreams
                [width, height, rate, nrFramesCaptured, nrFramesTotal, totalDuration] = mexDDGrab('getVideoInfo',i-1);
                video(i).width = width;
                video(i).height = height;
                video(i).rate = rate;
                video(i).nrFramesTotal = nrFramesTotal;
                video(i).totalDuration = totalDuration;
                video(i).frames = struct('cdata',cell(1,nrFramesCaptured),'colormap',cell(1,nrFramesCaptured));
                video(i).times = zeros(size(video(i).frames));
                video(i).skippedFrames = [];

                if (nrFramesTotal > 0 && any(frames > nrFramesTotal))
                    warning('mmread:general',['Frame(s) ' num2str(frames(frames>nrFramesTotal)) ' exceed the number of frames in the movie.']);
                end

                scanline = ceil(width*3/4)*4; % the scanline size must be a multiple of 4.

                for f=1:nrFramesCaptured
                    [data, time] = mexDDGrab('getVideoFrame',i-1,f-1);

                    if (numel(data) ~= scanline*height)
                        if (numel(data) > 3*width*height)
                            if (~warned)
                                warning('mmread:general','dimensions do not match data size. Guessing badly...');
                                warned = true;
                            end
                            scanline = width*3;
                            data = data(1:3*width*height);
                        else
                            if (f == 1)
                                error('dimensions do not match data size. Too little data.');
                            else
                                warning(['dimensions do not match data size. Too little data for ' num2str(f) 'th frame.']);
                                continue;
                            end
                        end
                    end

                    % if there is any extra scanline data, remove it
                    data = reshape(data,scanline,height);
                    data = data(1:3*width,:);

                    % the data ordering is wrong for matlab images, so permute it
                    tmp = permute(reshape(data, 3, width, height),[3 2 1]);
                    % the images are also upside down and colors were backwards.
                    video(i).frames(f).cdata = tmp(end:-1:1,:,3:-1:1);
                    video(i).times(f) = time;
                end

                framerate = (max(video(i).times)-min(video(i).times))/nrFramesCaptured;
                if framerate > 0
                    video(i).skippedFrames = any(diff(video(i).times)>framerate*1.8) & abs(mean(diff(video(i).times))-framerate)/framerate<0.05;
                end

                % if frames are specified then make sure that the order is the same
                if (~isempty(frames) && nrFramesCaptured > 0)
                    [uniqueFrames, dummy, frameOrder] = unique(frames);
                    if (length(uniqueFrames) > nrFramesCaptured)
                        warning('mmread:general','Not all frames specified were captured.  Returning what was captured, but order may be different than specified.');
                        remainingFrames = frames(frames<=uniqueFrames(nrFramesCaptured));
                        [dummy, dummy, frameOrder] = unique(remainingFrames);
                    end

                    video(i).frames = video(i).frames(frameOrder);
                    video(i).times = video(i).times(frameOrder);
                end
            end
        end

        % loop through getting all of the audio data from each stream
        for i=1:nrAudioStreams
            [nrChannels, rate, bits, nrFramesCaptured, nrFramesTotal, subtype, totalDuration] = mexDDGrab('getAudioInfo',i-1);
            audio(i).nrChannels = nrChannels;
            audio(i).rate = rate;
            audio(i).bits = bits;
            audio(i).nrFramesTotal = nrFramesTotal;
            audio(i).totalDuration = totalDuration;
            audio(i).frames = cell(1,nrFramesCaptured);
            for f=1:nrFramesCaptured
                [data, time] = mexDDGrab('getAudioFrame',i-1,f-1);
                audio(i).frames{f} = data;
                audio(i).times(f) = time;
            end
            % combine the data across frames
            d = double(cat(1,audio(i).frames{:}));

            % rescale the data so that it is between -1.0 and 1.0
            if (subtype==0)
                %PCM formated data...
                switch (bits)
                    case {4, 8}
                        d = (d-2^(bits-1))/2^(bits-1);
                    case {16, 24, 32}
                        d = d/2^(bits-1);
                end
            elseif (subtype==1)
                if (bits == 32)
                    %IEEE FLOAT formated data...
                    if (max(d) > 1 | min(d) < -1)
                        % there are two float formats one that is already -1 to 1
                        % and the there is between -2^15 to 2^15
                        d = d / 2^15;
                    end
                else
                    warning('Audio data format not recognized/supported, it probably is going to be useless.');
                end
            else
                warning('Audio data format not recognized/supported, it probably is going to be useless.');
            end

            % reshape the data so that it is nrChannels x Samples.  This should be the same output as wavread.
            audio(i).data = reshape(d,nrChannels,length(d)/nrChannels)';
        end

        mexDDGrab('cleanUp');
    catch
        err = lasterror;
   %     mexDDGrab('cleanUp');
        if strfind(err.message,'combination')
            disp('The ''No combination of intermediate filters could be found to make the connection'' error');
            disp('means that no appropriate codec could be found.  Mpg2 files seem to be the worst.  ');
            disp('Installing ffdshow (www.free-codecs.com/FFDShow_download.htm) often fixes this problem. ');
        end
        rethrow(err);
    end
end
