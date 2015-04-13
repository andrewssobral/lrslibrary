% ------------------------------------
% Copyright(c) <2009-2010> <Yin Zhang>
%   version 1.0, July 25, 2010
% ------------------------------------
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% -------------------------------
%   version 1.0, August 1, 2010
% -------------------------------
%
function [X,Y,S,Out] = lmafit_sms_v1(D,k,opts,beta)
%
% Low-rank matrix fitting for sparse matrix separation:
%                  LR + S ~= D
% where LR = XY is low-rank and S is sparse.  
% This code applies ADM to solving the model 
%         min  ||S||_1  s.t. XY + S = D.
%        X,Y,S
% Output:
%           X --- m x k matrix
%           Y --- k x n matrix
%           S --- m x n (sparse) matrix
%         Out --- structure for others
% Input:
%           D --- m x n data matrix
%           k --- initial rank estimate
%        opts --- structure for option with fields: 
%                 tol, maxit, print, ...
%        beta --- penalty parameter: a non-decreasing
%                 positive sequence of length >= 1               
%
[m,n] = size(D);
set_param(nargin);

% initialze
X = zeros(m,k);
Y = zeros(k,n);
S = zeros(m,n);
XY = S; L = S;

if iprint == 1; fprintf('Iteration:     '); end
if iprint == 2; RES = zeros(maxit,1); OBJ = RES; end

for iter = 1:maxit

    betav = beta(min(iter,nbeta));

    % min||XY - (D + L/beta - S)||_F^2 wrt X, then Y
    DL = D + (1/betav)*L;
    DLS = DL - S;
    if est_rank
         [X,R,E] = qr(DLS*Y',0);
    else [X,R  ] = qr(DLS*Y',0);
    end
    Y = X'*DLS;

    % save previous iterates
    if iter > 10 || iprint == 2; reso = resc;
        i1 = randi(m*n,5*m,1); XYo = XY(i1);
        i2 = randi(m*n,5*m,1);  So =  S(i2);
    end

    % shrink S: min_S||S - (D + L/beta - XY)||_F^2
    XY = X*Y;
    S = DL - XY;
    S = sign(S).*max(abs(S)-1/betav,0);
    
    % update multiplier
    Res = XY + S - D;
    L = L - betav*Res;

    % check stopping
    if iter > 10 || iprint == 2
        XYchg = norm(XYo-XY(i1))/norm(XYo);
        Schg  = norm(So - S(i2))/norm(So);
        res1 = norm(Res(i1))/norm(D(i1));
        res2 = norm(Res(i2))/norm(D(i2));
        resc = max(res1,res2);
        Reschg = abs(resc - reso)/max(1,reso);
        stop = max([XYchg Schg Reschg]) < tol;
        if stop; break; end
    end
    if iprint == 1; fprintf('\b\b\b\b\b%5i',iter); end
    if iprint == 2; 
        RES(iter) = res; OBJ(iter) = norm(S(:),1);
        fprintf('iter %4i: k %i XYchg %.2e Schg %.2e Res %.2e\n',...
            iter,k,XYchg,Schg,RES(iter));
    end
    
    % rank estimation
    if ~exist('rk','var'); rk = k; cnt = 0; icu = 0; end
    if est_rank && iter > 2 && k > 2; rank_estimator; end
    if rk < k;
        if iprint == -1
            fprintf('iter %i: est. rank %i --> %i\n',iter,k,rk);
        end
        k = max(rk,min_rank); 
        X = X(:,1:k); Y = Y(1:k,:); 
        est_rank = false;       
    end
        
end
X = X*sqrt(Dscale); Y = Y*sqrt(Dscale); S = S*Dscale;
if transposed; T = X; X = Y'; Y = T'; S = S'; end
    
if iprint == 1; fprintf('\n'); end
    Out.iter = iter;
    Out.rank = rk;
if iprint == 2;
    Out.res = RES(1:iter);
    Out.obj = OBJ(1:iter);
    Out.dual = DL - D;
end

% dammy statement
if 0 > 1; disp(R(1,1)+E(1,1)); end


%% nested function %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    function set_param(nargin)
        % process optional input
        tol = 1e-4;
        maxit = 500;
        iprint = 0;
        est_rank = 1;
        min_rank = min(1,k);
        if nargin >= 3
            if isfield(opts,'tol');       tol = opts.tol;    end
            if isfield(opts,'maxit');   maxit = opts.maxit;  end
            if isfield(opts,'print');  iprint = opts.print;  end
            if isfield(opts,'est_rank'); est_rank = opts.est_rank; end
            if isfield(opts,'min_rank'); min_rank = opts.min_rank; end
        end
        % transpose D if necessary
        transposed = false;
        if m >= 5*n && n >= 200
            D = D'; [m,n] = size(D);
            transposed = true;
        end
        Dscale = norm(D(:),inf);
        D = D / Dscale;
        if nargin < 2 || isempty(k); 
            k = round(0.5*min(m,n)); 
        end
        % set the beta sequence
        if nargin < 4 || isempty(beta);
            om = floor(log10(1/tol));
            nst = min(50,5*om);
            beta = 3/om*1.5.^(0:nst);
        else
            beta = beta*Dscale;
        end
        nbeta = length(beta);
        resc = inf;
         
    end %set_param

    %%%%%%%%%%%%%%%%%%%%%%%%
    function rank_estimator
        
        % parameters %%%%%%
        vbig = 0.95; % < 1
        vmed = 0.75; % 
        large = 3;   % >= 2
        max_cnt = 5; % > 1
        %%%%%%%%%%%%%%%%%%%
        
        mk1 = min_rank - 1;
        dR = abs(diag(R));
        wd = wdiff(dR);
        if mk1 > 0; wd(1:mk1) = 0; end
        [dmax,imaxd] = max(wd);
        if dmax > vbig;
            rk = imaxd; %figure(8);bar(wd);%pause;      
            return;
        end
        wr = wratio(dR);
        if mk1 > 0; wr(1:mk1) = 0; end
        [rmax,imaxr] = max(wr);
        ibig = find(wr > large);
        ibig2 = find(wr > .25*rmax);
        if numel(ibig) == 1 || numel(ibig2) == 1
            rk = imaxr; %figure(9);bar(wr);%pause;
            return;
        end
        wd2 = wdiff(wd).*(wd(1:end-1)>.1);
        if mk1 > 0; wd2(1:mk1)= 0; end
        [dmax2,imaxd2] = max(wd2);
        if dmax2 > vmed;
            rk = imaxd2; 
            return;
        end
        if wr(imaxd) > large;
            rk = imaxd;  
            return;
        end
        if imaxd == icu;
            cnt = cnt + 1;
        else
            icu = imaxd; cnt = 1; 
        end
        if cnt == max_cnt;
            rk = icu;
        end

    end %est rank

end %main

%%%%%%%%%%%%%%%%%%%%%%
function u = wdiff(v)
% input : v -- positive n-vector
% output: u -- weighted-diff (n-1)-vector
u = cumsum( v(end:-1:1) );
u = -diff(v) ./ u(end:-1:2);
end


%%%%%%%%%%%%%%%%%%%%%%
function u = wratio(v)
% input : v -- positive n-vector
% output: u -- weighted-ratio, (n-1)-vector
n = length(v);
u = v(1:end-1)./v(2:end);
u = (n-2)*u./(sum(u) - u);
end
