data = 'lobby.mat'; 
% 'mall.mat'
% 'yaleB02.mat'
% 'lobby.mat'
% 'hall.mat'
fprintf('--------------------------------\n')
fprintf(strcat(data, ' experiment', ' has started! \n'))
path = strcat('..\FW_SPCP\data\',data);
load(path); 
[m n] = size(D); 
D = D/norm(D,'fro');

%% parameter tuning
rho = 0.5;  % sampling ratio
Omega = rand(m,n)<=rho; % support of observation
obs = Omega.*D; % measurements

delta = 0.01;
% this is parameter to control noise level
% the smaller the noise, the smaller is delta

lambda_1 = delta*rho*norm(obs,'fro'); 
lambda_2 = delta*sqrt(rho)*norm(obs,'fro')/sqrt(max(m,n));

% compare or not
compare = 1;

par.M = obs; 
par.lambda_1 = lambda_1; par.lambda_2 =lambda_2;
par.iter = 1000; par.display = 1; par.rho = rho;
par.epsilon = 10^-3; % stopping criterion
par.method = 'exact';%'exact' or 'power'
par.Omega = Omega;
par.compare = compare; % make comparison or not
output = FW_T(par); % main function

% obtain the objective value returned from FW-T
L = output.L; S = output.S;
obj = 0.5*norm(Omega.*(L+S-obs),'fro')^2 + lambda_1*sum(svd(L))+...
    lambda_2*norm(vec(S),1);

if compare % make comparison with FISTA and ISTA
    par.obj = obj; 
    fprintf('the objective value returned by FW-T is %d \n', obj)
    tic
    fprintf('--------------------------------\n')
    fprintf('fista is invoked! \n')
    fista(par);
    toc
    tic
    fprintf('--------------------------------\n')
    fprintf('ista is invoked! \n')
    ista(par);
    toc
end
% plot objective values
figure();
plot(output.hist(1:end)); 
xlabel('iter.'); ylabel('obj. value');
title('obj. values produced by FW-T');

% original video
toplay = Omega.*D;
temp = size(toplay);
no_frames = temp(2); 
width = frameSize(1);
length = frameSize(2);
video_ori = zeros(width,length,no_frames); % width*length*frames
for i=1:no_frames
   video_ori(:,:,i) = reshape(toplay(:,i),width,length); 
end
% implay(video_ori/256,100);


% background term video
toplay = output.L;
no_frames = size(toplay,2);
width = frameSize(1);
length = frameSize(2);
video_L = zeros(width,length,no_frames); % width*length*frames
for i=1:no_frames
   video_L(:,:,i) = reshape(toplay(:,i),width,length); 
end

% sparse term video
toplay = output.S;
no_frames = size(toplay, 2);
width = frameSize(1);
length = frameSize(2);
video_S = zeros(width,length,no_frames); % width*length*frames
for i=1:no_frames
   video_S(:,:,i) = reshape(toplay(:,i),width,length); 
end
figure();
video_combine = [video_ori, video_L, abs(video_S)];
for i=1:no_frames    
    imagesc(video_combine(:,:,i));
    colormap('gray'); axis off;  title('FW-T');
    pause(0.01);
end



