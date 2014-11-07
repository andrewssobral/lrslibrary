function [eigvector, eigvalue, elapse] = LPP(W, options, data)
% LPP: Locality Preserving Projections
%
%       [eigvector, eigvalue] = LPP(W, options, data)
% 
%             Input:
%               data       - Data matrix. Each row vector of fea is a data point.
%               W       - Affinity matrix. You can either call "constructW"
%                         to construct the W, or construct it by yourself.
%               options - Struct value in Matlab. The fields in options
%                         that can be set:
%                           
%                         Please see LGE.m for other options.
%
%             Output:
%               eigvector - Each column is an embedding function, for a new
%                           data point (row vector) x,  y = x*eigvector
%                           will be the embedding result of x.
%               eigvalue  - The sorted eigvalue of LPP eigen-problem. 
%               elapse    - Time spent on different steps 
% 
%
%    Examples:
%
%       fea = rand(50,70);
%       options = [];
%       options.Metric = 'Euclidean';
%       options.NeighborMode = 'KNN';
%       options.k = 5;
%       options.WeightMode = 'HeatKernel';
%       options.t = 5;
%       W = constructW(fea,options);
%       options.PCARatio = 0.99
%       [eigvector, eigvalue] = LPP(W, options, fea);
%       Y = fea*eigvector;
%       
%       
%       fea = rand(50,70);
%       gnd = [ones(10,1);ones(15,1)*2;ones(10,1)*3;ones(15,1)*4];
%       options = [];
%       options.Metric = 'Euclidean';
%       options.NeighborMode = 'Supervised';
%       options.gnd = gnd;
%       options.bLDA = 1;
%       W = constructW(fea,options);      
%       options.PCARatio = 1;
%       [eigvector, eigvalue] = LPP(W, options, fea);
%       Y = fea*eigvector;
% 
% 
% Note: After applying some simple algebra, the smallest eigenvalue problem:
%				data^T*L*data = \lemda data^T*D*data
%      is equivalent to the largest eigenvalue problem:
%				data^T*W*data = \beta data^T*D*data
%		where L=D-W;  \lemda= 1 - \beta.
% Thus, the smallest eigenvalue problem can be transformed to a largest 
% eigenvalue problem. Such tricks are adopted in this code for the 
% consideration of calculation precision of Matlab.
% 
%
% See also constructW, LGE
%
%Reference:
%	Xiaofei He, and Partha Niyogi, "Locality Preserving Projections"
%	Advances in Neural Information Processing Systems 16 (NIPS 2003),
%	Vancouver, Canada, 2003.
%
%   Xiaofei He, Shuicheng Yan, Yuxiao Hu, Partha Niyogi, and Hong-Jiang
%   Zhang, "Face Recognition Using Laplacianfaces", IEEE PAMI, Vol. 27, No.
%   3, Mar. 2005. 
%
%   Deng Cai, Xiaofei He and Jiawei Han, "Document Clustering Using
%   Locality Preserving Indexing" IEEE TKDE, Dec. 2005.
%
%   Deng Cai, Xiaofei He and Jiawei Han, "Using Graph Model for Face Analysis",
%   Technical Report, UIUCDCS-R-2005-2636, UIUC, Sept. 2005
% 
%	Xiaofei He, "Locality Preserving Projections"
%	PhD's thesis, Computer Science Department, The University of Chicago,
%	2005.
%
%   version 2.1 --June/2007 
%   version 2.0 --May/2007 
%   version 1.1 --Feb/2006 
%   version 1.0 --April/2004 
%
%   Written by Deng Cai (dengcai2 AT cs.uiuc.edu)
%

bGlobal = 0;
if ~exist('data','var')
    bGlobal = 1;
    global data;
end

if (~exist('options','var'))
   options = [];
end

[nSmp,nFea] = size(data);
if size(W,1) ~= nSmp
    error('W and data mismatch!');
end


%==========================
% If data is too large, the following centering codes can be commented 
% options.keepMean = 1;
%==========================
if isfield(options,'keepMean') & options.keepMean
    ;
else
    if issparse(data)
        data = full(data);
    end
    sampleMean = mean(data);
    data = (data - repmat(sampleMean,nSmp,1));
end
%==========================




D = full(sum(W,2));


if ~isfield(options,'Regu') | ~options.Regu
    DToPowerHalf = D.^.5;
    D_mhalf = DToPowerHalf.^-1;

    if nSmp < 5000
        tmpD_mhalf = repmat(D_mhalf,1,nSmp);
        W = (tmpD_mhalf.*W).*tmpD_mhalf';
        clear tmpD_mhalf;
    else
        [i_idx,j_idx,v_idx] = find(W);
        v1_idx = zeros(size(v_idx));
        for i=1:length(v_idx)
            v1_idx(i) = v_idx(i)*D_mhalf(i_idx(i))*D_mhalf(j_idx(i));
        end
        W = sparse(i_idx,j_idx,v1_idx);
        clear i_idx j_idx v_idx v1_idx
    end
    W = max(W,W');
    
    data = repmat(DToPowerHalf,1,nFea).*data;
    [eigvector, eigvalue, elapse] = LGE(W, [], options, data);
else
    options.ReguAlpha = options.ReguAlpha*sum(D)/length(D);

    D = sparse(1:nSmp,1:nSmp,D,nSmp,nSmp);
    
    if bGlobal & isfield(options,'keepMean') & options.keepMean
        [eigvector, eigvalue, elapse] = LGE(W, D, options);
    else
        [eigvector, eigvalue, elapse] = LGE(W, D, options, data);
    end
end


eigIdx = find(eigvalue < 1e-3);
eigvalue (eigIdx) = [];
eigvector(:,eigIdx) = [];



