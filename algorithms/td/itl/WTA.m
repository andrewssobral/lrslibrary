function [T, C] = WTA(Xnew, R, Xold, overlap, type, C)
% Window-based tensor decomposition WTA 
%   
%   [T, C] = WTA((Xnew, R, Xold, overlap, type, C) compute window-based tensor
%   decomposition, according to  Xnew (Xold) the new (old) tensor window,
%   overlap is the number of overlapping tensors and C variance matrices
%   except for the time mode of previous tensor window and the specified
%   dimensions in vector R, type can be 'tucker' (default) or 'parafac'.
%   The input Xnew (Xold) is a tensor,sptensor, where the first mode is time.
%
%   The result returned in T is a tucker or kruskal tensor depending on
%   type and C is the cell array of new variance matrices.
%
%   Examples:
%
%----------------------------------------------------------------
% Copyright: 2006
%         Jimeng Sun, Christos Faloutsos
%         All rights reserved.
% Please address questions or comments to: jimeng@cs.cmu.edu
%----------------------------------------------------------------
if nargin<2
    error('not enough inputs');
end
N = ndims(Xnew);
if nargin<4
    overlap = 0;
end
if nargin<5
   type = 'tucker';
end

if nargin<6
    C = {};
    dv = size(Xnew);
    for n = 2:N
        C{n} = sparse(dv(n),dv(n));
    end
end
if strcmp(type, 'parafac') %for implementation convenience
    R = R*ones(N,1);
end
wsizenew = size(Xnew,1);
if exist('Xold')
    wsizeold = size(Xold,2);
end
if overlap > wsizenew 
    error('overlap size %d > Xnew size %d', overlap, wsizenew);
elseif exist('Xold') && overlap > wsizeold
    error('overlap size %d > Xold size %d', overlap, wsizeold);
elseif overlap < 0
    error('overlap size %d < 0', overlap);
end

%% Find the newly added and deleted tensors
if overlap == 0 
    Xadd = Xnew;
    if exist('Xold')
        Xdel = Xold;
    end
else
    Xadd = ['Xnew(', int2str(overlap+1),':end'];
    for n = 2:N
        Xadd = [Xadd,',:'];
    end
    Xadd = [Xadd,')'];
    Xadd = eval(Xadd);

    Xdel = ['Xold(1:', int2str(overlap)];
    for n = 2:N
        Xdel = [Xdel,',:'];
    end
    Xdel = [Xdel,')'];
    Xdel = eval(Xdel);
end
%% Update covariance matrices C and compute U for initialization
U = {};
for n = 2:N %every mode except time (1st mode)
    switch class(Xnew)
        case {'tensor'} %dense
            XMadd = double(tenmat(Xadd,n));
        case {'sptensor'} %sparse
            XMadd = double(sptenmat(Xadd,n));
    end
    C{n} = C{n} + XMadd*XMadd';
    if exist('Xold')
        switch class(Xold)
            case {'tensor'} %dense
                XMdel = double(tenmat(Xdel,n));
            case {'sptensor'} %sparse
                XMdel = double(sptenmat(Xdel,n));
        end
        C{n} = C{n} - XMdel*XMdel';    
    end

    opts.disp = 0;
    opts.issym = 1;
    warning off;
    [U{n}, D] = eigs(C{n},R(n),'LM',opts);
    warning on;
end

%% tensor decomposition
if strcmp(type,'parafac')
    T = parafac_als(Xnew, R(1), struct('dimorder',1:N,'init',{U}));
elseif strcmp(type,'tucker')
    T = tucker_als(Xnew, R, struct('dimorder',1:N,'init',{U}));
end
