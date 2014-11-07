function [img2warp,tau] = regMGNC(img1,img2,tau,numLevel,paraTurkey0)

img1 = double(img1);
img2 = double(img2);

if nargin < 3 || isempty(tau)
    tau = zeros(6,1);
end
if nargin < 4 || isempty(numLevel)
    numLevel = 3;
end
if nargin < 5 || isempty(paraTurkey0)
    % that means minimum C depends on the online computing
    paraTurkey0 = eps;
end

for l = numLevel-1:-1:0
    %% construct pyramid
    I1 = imresize(img1,0.5^l,'Antialiasing',true);
    I2 = imresize(img2,0.5^l,'Antialiasing',true);
    
    % initialize
    if l == numLevel-1
        weight = ones(size(I1));
        tau(end-1:end) = 0.5^(l-1)*tau(end-1:end); % pass from the initial value
        C = max(abs(img1(:))); % the coarsest level, C starts from a big number
    else
        weight = imresize(weight,size(I1));
        tau(end-1:end) = 2*tau(end-1:end); % pass from previous level
        C = C*2;
    end
    
    while 1
       %% gradually nonconvexity
       
        % reweighted least square
        for iter = 1:50
            tau_old = tau;
            [img2warp,tau,residue] = regImg(I1,I2,tau_old,weight,1);
            weight = (influence(abs(residue),C)+eps)./(abs(residue)+eps); 
            weight = weight/max(weight(:));     
            % check termination condition
            if max(abs((tau_old-tau)./(tau+eps))) < 0.01;
                break
            end
        end  
        
        % estimate paraTurkey only
        paraTurkey = 4.7*1.48* median(abs( residue(:) - median(residue(:)) ));
        
        % check C
        if C >= max( paraTurkey,paraTurkey0 )
            C = C/2;
        else
            break;
        end 
     
    end
   
end

end

function y = influence(x,C)
% Influence function from Tukey's biweight estimator
y = zeros(size(x));
idx = abs(x)<C;
y(idx) = x(idx).*(C^2-x(idx).^2).^2;
end