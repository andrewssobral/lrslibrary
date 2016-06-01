function [X,Y,Out] = lmafit_mc_adp(m,n,k,Known,data,opts)
%
% Solver for matrix completion:
%
%        U_ij ~= A_ij,  for i,j in Known
%
% where U =  X*Y and A_ij's are given in an known index set.
%
% Output:
%           X --- m x k matrix
%           Y --- k x n matrix
%         Out --- output information
% Input:
%        m, n --- matrix sizes
%           k --- rank estimate
%        data --- values of known elements in a 1D row vector
%       Known --- positions of known elements in a 1D row vector
%                 assuming matrices are arranged column-wise
%             or  Known is a structure with two fields [Ik, Jk]
%        opts --- option structure with fields: opt, maxit
%
% Copyright(c) 2009 Yin Zhang
%
%   modified by Zaiwen Wen, 12/17/2009
%

L = length(data);

% set parameters
tol = 1.25e-4;
maxit = 500;
iprint = 1;
Zfull = (L/(m*n) > 0.2 ) || k > .02*min(m,n) || m*n < 5e5;
DoQR = true;
est_rank = 1;
rank_max =  floor(0.1*min(m,n));
rank_min =  1;
rk_inc = 5;
rk_jump = 10;
init = 0;
save_res = 0;

if isfield(opts,'tol');         tol     = opts.tol;        end
if isfield(opts,'maxit');       maxit   = opts.maxit;      end
if isfield(opts,'print');       iprint  = opts.print;      end
if isfield(opts,'Zfull');       Zfull   = opts.Zfull;      end
if isfield(opts,'DoQR');        DoQR    = opts.DoQR;       end
if isfield(opts,'est_rank');    est_rank= opts.est_rank;   end
if isfield(opts,'rank_max');    rank_max= opts.rank_max;   end
if isfield(opts,'rank_min');    rank_min= opts.rank_min;   end
if isfield(opts,'rk_inc');      rk_inc  = opts.rk_inc;     end 
if isfield(opts,'rk_jump');     rk_jump = opts.rk_jump;    end 
if isfield(opts,'init');        init    = opts.init;       end
if isfield(opts,'save_res');    save_res= opts.save_res;   end



reschg_tol = 0.5*tol; 
rk = k;  

% if est_rank == 1; rank_max = min(rank_max, k); end; 
    
linopts.SYM = true; linopts.POSDEF = true;
datanrm = max(1,norm(data));    

objv = zeros(maxit,1);  RR = ones(maxit,1);
if iprint == 1; fprintf('Iteration:     '); end
if iprint == 2;
    fprintf('\nLMafit_mc: Zfull = %i, DoQR = %i\n',Zfull,DoQR);
end

% initialize: make sure the correctness of the index set and data
data(data==0) = eps;    data_tran = false; 

Zfull = 0; % to be studied

if Zfull %Z is full
    Z = zeros(m,n); 
    Z(Known) = data;   
else %Z = S + XY, initialize the storage of S
    if isnumeric(Known);       [Ik,Jk] = ind2sub([m n],Known);
    elseif isstruct(Known)     Ik = Known.Ik; Jk = Known.Jk;     end
    %make sure the order of Ik, Jk and data are correctly as in a sparse matrix
    S = sparse(Ik, Jk, data, m, n); 
    [Ik, Jk, data] = find(S);  
    data = data';   
end

if m>n; 
    tmp = m; m = n; n = tmp; data_tran = true; 
    if Zfull;  Z = Z';  Known = find(Z);        data = Z(Known);
    else       S = S';  [Ik, Jk, data] = find(S);  data = data';    end
end

if init == 0
    X = zeros(m,k);   Y = eye(k,n);   Res = data;   res = datanrm; 
    %X = zeros(m,k);   Y = rand(k,n);   Res = data;   res = datanrm; 
else
    X = opts.X;  Y = opts.Y;    opts.X = []; opts.Y = [];
    if Zfull 
        Z = X*Y;  Res = data - Z(Known);    Z(Known) = data;
    else %Z = S + XY, initialize the storage of S
        Res = data - partXY(X',Y,Ik,Jk,L);  updateSval(S, Res, L);
    end
    res = norm(Res);
end

% parameters for alf
alf = 0;  increment = 1; %init_inc();
itr_rank = 0; minitr_reduce_rank = 5;    maxitr_reduce_rank = 50;  

Out.newRank = [];
Out.stopCrtr =[];
 t0 = tic;

 
 % main iteration
for iter = 1:maxit
    
    itr_rank = itr_rank + 1; 
    Xo = X; Yo = Y; Res0 = Res; res0 = res; alf0x = alf;
    
    if Zfull
        Zo = Z; X = Z*Y';
        if est_rank == 1
            [X,R,E] = qr(X,0);  Y = X'*Z;
        elseif DoQR
            [X,R  ] = qr(X,0);  Y = X'*Z;
        else
            Xt = X'; Y = linsolve(Xt*X,Xt*Z,linopts);
        end
        Z = X*Y;     Res = data - Z(Known); 
    else % Z=S+XY for sparse S(Known)=data-XY(Known)
        Yt = Y';     X = S*Yt + X*(Y*Yt);    % Z*Y'
        if est_rank == 1
            [X,R,E] = qr(X,0);  Xt = X';
            Y = Xt*S + (Xt*Xo)*Y;       % X'*Z
        elseif DoQR
            [X,R  ] = qr(X,0);  Xt = X';
            Y = Xt*S + (Xt*Xo)*Y;
        else Xt = X';
            Y = Xt*S + (Xt*Xo)*Y;
            Y = linsolve(Xt*X,Y,linopts);
        end
        Res = data - partXY(Xt,Y,Ik,Jk,L);
    end %Zfull

    res = norm(Res);  relres = res/datanrm;  ratio = res/res0;
    reschg = abs(1-res/res0); RR(iter) = ratio; 
    

    % adjust alf
    if ratio >= 1 
        increment = max(0.1*alf, 0.1*increment);
        %increment = max(0.1*alf, 0.5*increment);
        X = Xo; Y = Yo; Res = Res0; res = res0; relres = res/datanrm;
        alf = 0;    if Zfull; Z = Zo;end
    elseif ratio > 0.7;
        increment = max(increment, 0.25*alf);   
        alf = alf + increment;
    end

    obj_val = res^2/L;
    % printout
    if iprint == 1; fprintf('\b\b\b\b\b%5i',iter); end
    if iprint == 2;
        fprintf('it: %5i obj. %3e rel. %3.1e r. %4.4f chg: %3.1e alf: %3.1e inc: %3.1e\n',...
            iter,obj_val, relres,ratio,reschg,alf0x,increment);
    end
    
    
    objv(iter) = obj_val;  % Mean square error
        
    
    % check stopping
%     if ((reschg < reschg_tol && ...
%             itr_rank > minitr_reduce_rank) || relres < tol); Out.stopCrtr=[Out.stopCrtr;relres]; break; end
    
    if obj_val <= tol,
        break;
    end;
    
    
    
    %if ((est_rank ==1||est_rank==0) && ((reschg < reschg_tol && ...
    %        itr_rank > minitr_reduce_rank) || relres < tol) ) ...
    %   || ((est_rank == 2)&&( relres < tol || (k==rank_max && itr_rank >= maxit)))
    %        break; 
    %end


    % update Z or S with the most recent alf
    if Zfull; Z(Known) = data + alf*Res; else updateSval(S, (alf+1)*Res, L); end
end %iter


if iprint == 1; fprintf('\n'); end
if data_tran; tX = X; X = Y'; Y = tX'; end
Out.obj = objv(1:iter);
Out.RR = RR(1:iter);
Out.iter = iter;
Out.rank = rk;
Out.relres = relres;
Out.reschg = reschg;
Out.time = toc(t0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     function init_inc()
%         dr = L/(m+n-k)/k;  increment = .5*log10(m*n) - log10(rk);
%         if min(m,n)/k > 1000; increment = increment + .25*exp(dr-1)/dr; end
%     end

%     function rank_estimator_adaptive()
%         if est_rank == 1
%             dR = abs(diag(R));       drops = dR(1:end-1)./dR(2:end);
%             [dmx,imx] = max(drops);  rel_drp = (k-1)*dmx/(sum(drops)-dmx);
%             if (rel_drp > rk_jump && itr_rank > minitr_reduce_rank) ...
%                    || itr_rank > maxitr_reduce_rank; %bar(drops); pause;
%                 rk = max([imx, floor(0.1*k), rank_min]); 
%                 X = X(:,1:rk); Y = Y(1:rk,:);
%                 if Zfull %Z is full
%                     Z = X*Y;  Res = data - Z(Known);   
%                     Z(Known) = data + alf*Res; 
%                 else %Z = S + XY
%                     Res = data - partXY(X',Y,Ik,Jk,L);
%                     updateSval(S, (alf+1)*Res, L); 
%                 end
%                 res = norm(Res);  est_rank = 0; itr_rank = 0; 
%                 if iprint >= 2 %&& ch_rank>=1
%                     fprintf('it: %5i, rel_drp: %3.2e, est_rank: %d,  Rank estimate changed from %i to %i\n',...
%                         iter, rel_drp, est_rank, k,rk);
%                 end
%             end
%         elseif est_rank == 2 && reschg < 10*tol && rk < rank_max && itr_rank > 1 
%         %elseif est_rank == 2 && reschg < tol/10  && itr_rank > 1 
%             
%             Out.newRank=[Out.newRank; iter]; % count the iterations for a rank change
%             Out.stopCrtr=[Out.stopCrtr; relres];
%             
%             if rk < 50; rinc = rk_inc; else  rinc = 2*rk_inc; end
%             rk = min(rk + rinc, rank_max); rkr_id = true;
%             if rk > k;
%                 if save_res == 1
%                     save(strcat('LM-Med-r', num2str(rk),'max', num2str(rank_max),  '.mat'), 'X', 'Y');
%                 end
%                 X = [X, zeros(m, rk-k)]; Y = [Y; zeros(rk-k,n)]; itr_rank = 0;
%                 if iprint >= 2 %&& ch_rank>=1
%                     fprintf('it: %5i, reschg: %3.2e, Rank estimate changed from %i to %i\n',...
%                         iter,reschg, k,rk);
%                 end
%             end
%         end
%     end %rank

end %main
