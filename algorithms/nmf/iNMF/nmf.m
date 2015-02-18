% Hoer's code (P. O. Hoyer. Non-negative Matrix Factorization with sparseness constraints. 
% Journal of Machine Learning Research  5:1457-1469, 2004.) is modified 

function [W, H, objhistory, objhistory_v] = nmf( V, rdim, showflag, ttt )
%% INPUT
% V: data matrix
% rdim : matrix factorization rank
% showflag: if it is 1. than it plots the change of reconstruction error
% for each iteration
% ttt: maximum number of iterations allowed

%% OUTPUT
% W: Mixing matrix (V=WH)
% H: Encoding matrix (V=WH)
% objhistory: change of total reconstruction error
% objhistory_v: reconstruction error for each individual sample




% Check that we have non-negative data
if min(V(:))<0, error('Negative values in data!'); end
 
% Dimensions
vdim = size(V,1);
samples = size(V,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% initialization %%%
fname2=['initials/Initials' num2str(vdim) 'x' num2str(rdim) '.txt'];
fid=fopen(fname2,'r');
if (fid==-1)
    ssbinitial(vdim,rdim);
end
fidW2 = fopen(fname2,'r');
W = fscanf(fidW2,'%f', [vdim inf]);
fclose(fidW2);

clear fname2

fname2=['initials/Initials' num2str(rdim) 'x' num2str(samples) '.txt'];
fid=fopen(fname2,'r');
if (fid==-1)
    ssbinitial(rdim,samples);
end
fidW2 = fopen(fname2,'r');
H = fscanf(fidW2,'%f', [rdim inf]);
fclose(fidW2);

clear fname2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (showflag==1)
    objhistory = ((sum(sum((V-W*H).^2))))/(vdim*samples);
    figure(2); clf;
    drawnow;
end


% Start iteration
iter = 0;
while 1,
    
    if (iter==ttt)
        break
    end
    iter = iter+1;    
    % Compute new W and H (Lee and Seung; NIPS*2000)
    H = H.*(W'*V)./(W'*W*H + 1e-9);
    W = W.*(V*H')./(W*H*H' + 1e-9);
    
%     norms = sqrt(sum(W.^2));
%     H = H.*(norms'*ones(1,samples));
%     W = W./(ones(vdim,1)*norms);

    newobj = ((sum(sum((V-W*H).^2))))/(vdim*samples);
    newobj_v = ((sum(sum((V(:,end)-W*H(:,end)).^2))))/(vdim);
    if iter==1
        objhistory = [newobj];
        objhistory_v = [newobj_v];
    else
        objhistory = [objhistory newobj];
        objhistory_v = [objhistory_v newobj_v];
    end
    
    if(showflag==1)
        if (iter>1)
            plot(objhistory(2:end)); 
        end
        drawnow;
    end   
end