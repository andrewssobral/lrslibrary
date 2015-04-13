function [ bgU, status, OPTS ] = bgtraining( videopath, OPTIONS, max_cycles, subsampling, training_size, FPS_ONLY , TRAIN_FRAME)
%BGTRAINING Summary of this function goes here
%   Detailed explanation goes here
thresh = 0.1;

FILE_EXT        = '.bmp';
OPTS            = struct(); % initial a empty struct for OPTS
status.init     = 0;   % status of grasta at each iteration
U_hat           = zeros(1);

%
Dir_lists = dir(videopath);
Video_Length = min(training_size,length(Dir_lists));

if TRAIN_FRAME == 1, % Random selection
    p = randperm(length(Dir_lists)); % length(Dir_lists)-1500 for hall
    frames_idx = p(1:Video_Length);
else                % Select all frames
    frames_idx = 1:Video_Length;
end

if ~FPS_ONLY,
    h_fg = subplot(2,2,1);set(gca,'nextplot','replacechildren');title('GRASTA-Foreground');
    h_training_bg = subplot(2,2,2);set(gca,'nextplot','replacechildren');title('Training Background');
    h_img = subplot(2,2,3);set(gca,'nextplot','replacechildren');title('Training video');
end

for outiter = 1:max_cycles,
    video_order = randperm(Video_Length);  %  1:Video_Length; % 
    t_start = tic;
    frame_count  = 0;
    for i=1:Video_Length,
        if Dir_lists(frames_idx(video_order(i))).isdir,
            continue;
        end
        if isempty(strfind(Dir_lists(frames_idx(video_order(i))).name,FILE_EXT)) ,
            continue;
        end
        
        frame_count = frame_count+1;
        fname = [videopath Dir_lists(frames_idx(video_order(i))).name];
        
        % prepare the image
        I = imread(fname);
        I = double(rgb2gray(I));
        
        if frame_count==1,
            [rows,cols]     = size(I);
            VW_ROWS         = rows;
            VW_COLS         = cols; %ceil(cols * VW_RATIO);
            DIM             = VW_ROWS * VW_COLS;
            OPTIONS.DIM_M   = DIM; % video ambient dimension
        end
        
        I = I/max(max(I));
                
        % random subsampling the frame I
        M = round(subsampling * DIM);
        p = randperm(DIM);
        idx = p(1:M)';
        
        I_Omega = I(idx);       
        
        [U_hat, status, OPTS] = grasta_stream(I_Omega, idx, U_hat, status, OPTIONS, OPTS);
        
        if mod(i,1) == 0 && ~FPS_ONLY,
            bg_img = reshape(U_hat * status.w * status.SCALE, VW_ROWS,VW_COLS);
            axes(h_training_bg); imagesc(bg_img);colormap gray;axis off;axis ij ;
            
            o_img = reshape( I ,VW_ROWS,VW_COLS );
            axes(h_img); imagesc(o_img);colormap gray;axis off;axis ij ;
            
            s_hat = I(:) - U_hat * status.w * status.SCALE;
            s_hat = fg_thresholding(s_hat,thresh);

            s_img = reshape(s_hat,VW_ROWS,VW_COLS);
            axes(h_fg); imagesc(s_img);colormap gray;axis off;axis ij ;
            
        end
    end
    t_end = toc(t_start);
    fprintf('Training %d/%d: %.2f seconds, %.2f fps, grasta_t %.2e \n',...
        outiter, max_cycles,t_end, frame_count/t_end,status.grasta_t);
end

bgU = U_hat;
end

