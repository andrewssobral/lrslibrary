%% double process_video(int, int, string, string);
% method_id - string
% algorithm_id - string
% inFile - string
% outFile - string
%
% demos:
% process_video('RPCA', 'FPCP', 'dataset/demo.avi', 'output/demo_out.avi');
% process_video('RPCA', 'FPCP', 'dataset/demo.avi', 'output/demo_out.mat');
% process_video('RPCA', 'FPCP', 'dataset/demo.mat', 'output/demo_out.mat');
% process_video('RPCA', 'FPCP', 'dataset/demo.mat', 'output/demo_out.avi');
%
% unix: 
% ./matlab -nojvm -nodisplay -nosplash -r "process_video('RPCA', 'FPCP', 'dataset/demo.avi', 'output/demo_out.avi');exit;"
% ./matlab -nojvm -nodisplay -nosplash -r "process_video('RPCA', 'FPCP', 'dataset/demo.mat', 'output/demo_out.avi');exit;"
%
% For debug:
% load('output/demo_out.mat');
% showResultsInfo(info);
%
% method_id='RPCA'; algorithm_id='FPCP'; inFile='dataset/demo.avi'; outFile='output/demo_out.avi';
% inFile = 'dataset/ChangeDetection2012/badminton_out.avi';
function [stats] = process_video(method_id, algorithm_id, inFile, outFile)
timerVal = tic;
displog(['Loading ' inFile]);
video = load_video_file(inFile);

%%% Matrix-based methods
% i.e: process_video('RPCA', 'FPCP', 'dataset/demo.avi', 'output/demo_FPCP.avi');
if(strcmp(method_id,'RPCA') || strcmp(method_id,'ST') || strcmp(method_id,'MC') ...
|| strcmp(method_id,'LRR') || strcmp(method_id,'TTD') || strcmp(method_id,'NMF'))
  M = im2double(convert_video_to_2d(video));
  params.rows = video.height;
  params.cols = video.width;
	results = run_algorithm(method_id, algorithm_id, M, params);
  movobj = convert_2dresults2mov([],results.L,results.S,results.O,video);
end

%%% Tensor-based methods
% i.e: process_video('TD', 'HOSVD', 'dataset/demo.avi', 'output/demo_HOSVD.avi');
if(strcmp(method_id,'TD') || strcmp(method_id,'NTF'))
  A = im2double(convert_video_to_3d(video));
  T = tensor(A);
	results = run_algorithm(method_id, algorithm_id, T, []);
  movobj = convert_3dresults2mov([],results.L,results.S,results.O,size(T,3));
end

displog('Saving results...');
save_results(movobj,outFile);

displog('Process finished!');
displog(['CPU time: ' num2str(results.cputime)]);

stats.cputime = results.cputime; % Elapsed time for decomposition
stats.totaltime = toc(timerVal); % Total elapsed time
end
