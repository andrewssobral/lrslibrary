%% This function is to seperate the given video frames into foreground/background 
% in realtime by our proposed algorithm - GRASTA, and return the corresponding 
% seperated results of ground-truths as fgs_cmp.
%
% Author: Jun He, Laura Balzano, Arthur Szlam
% Date:  Oct. 30, 2011
%%

function [fgs_cmp,bgs_cmp,vInfo] = bgfg_seperation_grasta( U_hat, videopath, subsampling , status, OPTIONS, OPTS,...
                                       thresh,FPS_ONLY, frame_names, MAX_FRAME)
                                   
FILE_EXT        = '.bmp';
frame_count     = 0;

%
Dir_lists = dir(videopath);
Video_Length = length(Dir_lists);

if ~FPS_ONLY,    
    figure;
    h_fg = subplot(2,2,1);set(gca,'nextplot','replacechildren');title('Foreground');
    h_fg_bw = subplot(2,2,2);set(gca,'nextplot','replacechildren');title('Thresholded-Foreground');
    h_bg = subplot(2,2,3);set(gca,'nextplot','replacechildren');title('Background');
    h_img = subplot(2,2,4);set(gca,'nextplot','replacechildren');title('Video');
end

t_start = tic;
for i=1:Video_Length,
    if Dir_lists(i).isdir,
        continue;
    end
    if isempty(strfind(Dir_lists(i).name,FILE_EXT)) ,
        continue;
    end
    
    frame_count = frame_count+1;
    fname = [videopath Dir_lists(i).name];
    
    % prepare the image
    I = imread(fname);
    I = double(rgb2gray(I));
    
    if frame_count==1,
        [rows,cols]     = size(I);
        VW_ROWS         = rows;
        VW_COLS         = cols; %ceil(cols * VW_RATIO);
        DIM             = VW_ROWS * VW_COLS;
        OPTIONS.DIM_M   = DIM; % video ambient dimension
        fgs_cmp         = zeros(DIM, length(frame_names));
        bgs_cmp         = zeros(DIM, length(frame_names));
        
        vInfo.rows = rows;
        vInfo.cols = cols;
    end
    
    I = I/max(max(I));
    
    % random subsampling the frame I
    if subsampling < 1,
        M = round(subsampling * DIM);
        p = randperm(DIM);
        idx = p(1:M)';
    else
        idx = (1:DIM)';
    end
    
    I_Omega = I(idx);
    
    % tracking the background
    [U_hat, status, OPTS] = grasta_stream(I_Omega, idx, U_hat, status, OPTIONS, OPTS);  
    
    bg_img = reshape(U_hat * status.w * status.SCALE, VW_ROWS,VW_COLS);
    o_img = reshape( I ,VW_ROWS,VW_COLS );   
    
    noise_thresh = 1 * min(abs(I(:)));
    s_hat = I(:) - U_hat * status.w * status.SCALE;
    s_hat(abs(s_hat) < noise_thresh) = 0;
    s_img = reshape(s_hat,VW_ROWS,VW_COLS);
    
    s_hat = fg_thresholding(s_hat,thresh);
    s_img_bw = reshape(s_hat,VW_ROWS,VW_COLS);
    
    for jj=1:length(frame_names),
        if strcmpi(Dir_lists(i).name, frame_names{jj}),
            fgs_cmp(:,jj) = s_img(:);
            bgs_cmp(:,jj) = bg_img(:);
            break;
        end
    end
    
    if ~FPS_ONLY,
        axes(h_bg); imagesc(bg_img);colormap gray;axis off;axis ij ;
        axes(h_img); imagesc(o_img);colormap gray;axis off;axis ij ;
        axes(h_fg); imagesc(s_img);colormap gray;axis off;axis ij ;        
        axes(h_fg_bw); imagesc(s_img_bw);colormap gray;axis off;axis ij ;        
    end
    
    if frame_count >= MAX_FRAME && MAX_FRAME~=-1,
        break;
    end
end
t_end = toc(t_start);
fprintf('Tracking %d frames with %.1f%% information costs %.2f seconds, %.2f fps\n',...
    frame_count, 100*subsampling,t_end, frame_count/t_end);
end

