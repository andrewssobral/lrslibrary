function [eigvector, eigvalue, elapse] = LDA(gnd,options,data)
% LDA: Linear Discriminant Analysis 
%
%       [eigvector, eigvalue] = LDA(gnd, options, data)
% 
%             Input:
%               data  - Data matrix. Each row vector of fea is a data point.
%               gnd   - Colunm vector of the label information for each
%                       data point. 
%               options - Struct value in Matlab. The fields in options
%                         that can be set:
%
%                            Regu  -  1: regularized solution, 
%                                        a* = argmax (a'X'WXa)/(a'X'Xa+ReguAlpha*I) 
%                                     0: solve the sinularity problem by SVD 
%                                     Default: 0 
%
%                         ReguAlpha -  The regularization parameter. Valid
%                                      when Regu==1. Default value is 0.1. 
%
%                            ReguType  -  'Ridge': Tikhonov regularization
%                                         'Custom': User provided
%                                                   regularization matrix
%                                          Default: 'Ridge' 
%                        regularizerR  -   (nFea x nFea) regularization
%                                          matrix which should be provided
%                                          if ReguType is 'Custom'. nFea is
%                                          the feature number of data
%                                          matrix
%                        Fisherface     -  1: Fisherface approach
%                                             PCARatio = nSmp - nClass
%                                          Default: 0
%
%                            PCARatio     -  The percentage of principal
%                                            component kept in the PCA
%                                            step. The percentage is
%                                            calculated based on the
%                                            eigenvalue. Default is 1
%                                            (100%, all the non-zero
%                                            eigenvalues will be kept.
%                                            If PCARatio > 1, the PCA step
%                                            will keep exactly PCARatio principle
%                                            components (does not exceed the
%                                            exact number of non-zero components).  
%                           
%
%             Output:
%               eigvector - Each column is an embedding function, for a new
%                           data point (row vector) x,  y = x*eigvector
%                           will be the embedding result of x.
%               eigvalue  - The sorted eigvalue of LDA eigen-problem. 
%               elapse    - Time spent on different steps 
%
%    Examples:
%       
%       fea = rand(50,70);
%       gnd = [ones(10,1);ones(15,1)*2;ones(10,1)*3;ones(15,1)*4];
%       options = [];
%       options.Fisherface = 1;
%       [eigvector, eigvalue] = LDA(gnd, options, fea);
%       Y = fea*eigvector;
% 
%
% See also LPP, constructW, LGE
%
%
%
%Reference:
%
%   P. N. Belhumeur, J. P. Hespanha, and D. J. Kriegman, “Eigenfaces
%   vs. fisherfaces: recognition using class specific linear
%   projection,” IEEE Transactions on Pattern Analysis and Machine
%   Intelligence, vol. 19, no. 7, pp. 711-720, July 1997.  
%
%   Deng Cai, Xiaofei He, Yuxiao Hu, Jiawei Han, and Thomas Huang, 
%   "Learning a Spatially Smooth Subspace for Face Recognition", CVPR'2007
%
%   Deng Cai, Xiaofei He, Jiawei Han, "SRDA: An Efficient Algorithm for
%   Large Scale Discriminant Analysis", IEEE Transactions on Knowledge and
%   Data Engineering, 2007.
%
%   version 2.1 --June/2007 
%   version 2.0 --May/2007 
%   version 1.1 --Feb/2006 
%   version 1.0 --April/2004 
%
%   Written by Deng Cai (dengcai2 AT cs.uiuc.edu)
%                                                   

if ~exist('data','var')
    global data;
end

if (~exist('options','var'))
   options = [];
end

if ~isfield(options,'Regu') | ~options.Regu
    bPCA = 1;
    if ~isfield(options,'PCARatio')
        options.PCARatio = 1;
    end
else
    bPCA = 0;
    if ~isfield(options,'ReguType')
        options.ReguType = 'Ridge';
    end
    if ~isfield(options,'ReguAlpha')
        options.ReguAlpha = 0.1;
    end
end

tmp_T = cputime;

% ====== Initialization
[nSmp,nFea] = size(data);
if length(gnd) ~= nSmp
    error('gnd and data mismatch!');
end

classLabel = unique(gnd);
nClass = length(classLabel);
Dim = nClass - 1;

if bPCA & isfield(options,'Fisherface') & options.Fisherface
    options.PCARatio = nSmp - nClass;
end


if issparse(data)
    data = full(data);
end
sampleMean = mean(data,1);
data = (data - repmat(sampleMean,nSmp,1));



bChol = 0;
if bPCA & (nSmp > nFea+1) & (options.PCARatio >= 1)
    DPrime = data'*data;
    DPrime = max(DPrime,DPrime');
    [R,p] = chol(DPrime);
    
    if p == 0
        bPCA = 0;
        bChol = 1;
    end
end


%======================================
% SVD
%======================================
if bPCA    
    if nSmp > nFea
        ddata = data'*data;
        ddata = max(ddata,ddata');

        [eigvector_PCA, eigvalue_PCA] = eig(ddata);
        eigvalue_PCA = diag(eigvalue_PCA);
        clear ddata;

        maxEigValue = max(abs(eigvalue_PCA));
        eigIdx = find(eigvalue_PCA/maxEigValue < 1e-12);
        eigvalue_PCA(eigIdx) = [];
        eigvector_PCA(:,eigIdx) = [];

        [junk, index] = sort(-eigvalue_PCA);
        eigvalue_PCA = eigvalue_PCA(index);
        eigvector_PCA = eigvector_PCA(:, index);
        
        %=======================================
        if options.PCARatio > 1
            idx = options.PCARatio;
            if idx < length(eigvalue_PCA)
                eigvalue_PCA = eigvalue_PCA(1:idx);
                eigvector_PCA = eigvector_PCA(:,1:idx);
            end
        elseif options.PCARatio < 1
            sumEig = sum(eigvalue_PCA);
            sumEig = sumEig*options.PCARatio;
            sumNow = 0;
            for idx = 1:length(eigvalue_PCA)
                sumNow = sumNow + eigvalue_PCA(idx);
                if sumNow >= sumEig
                    break;
                end
            end
            eigvalue_PCA = eigvalue_PCA(1:idx);
            eigvector_PCA = eigvector_PCA(:,1:idx);
        end
        %=======================================
        
        eigvalue_PCA = eigvalue_PCA.^-.5;
        data = (data*eigvector_PCA).*repmat(eigvalue_PCA',nSmp,1);
    else
        ddata = data*data';
        ddata = max(ddata,ddata');

        [eigvector, eigvalue_PCA] = eig(ddata);
        eigvalue_PCA = diag(eigvalue_PCA);
        clear ddata;

        maxEigValue = max(eigvalue_PCA);
        eigIdx = find(eigvalue_PCA/maxEigValue < 1e-12);
        eigvalue_PCA(eigIdx) = [];
        eigvector(:,eigIdx) = [];

        [junk, index] = sort(-eigvalue_PCA);
        eigvalue_PCA = eigvalue_PCA(index);
        eigvector = eigvector(:, index);
        
        %=======================================
        if options.PCARatio > 1
            idx = options.PCARatio;
            if idx < length(eigvalue_PCA)
                eigvalue_PCA = eigvalue_PCA(1:idx);
                eigvector = eigvector(:,1:idx);
            end
        elseif options.PCARatio < 1
            sumEig = sum(eigvalue_PCA);
            sumEig = sumEig*options.PCARatio;
            sumNow = 0;
            for idx = 1:length(eigvalue_PCA)
                sumNow = sumNow + eigvalue_PCA(idx);
                if sumNow >= sumEig
                    break;
                end
            end
            eigvalue_PCA = eigvalue_PCA(1:idx);
            eigvector = eigvector(:,1:idx);
        end
        %=======================================
        

        eigvalue_PCA = eigvalue_PCA.^-.5;
        eigvector_PCA = (data'*eigvector).*repmat(eigvalue_PCA',nFea,1);

        data = eigvector;
        clear eigvector;
    end
else
    if ~bChol
        DPrime = data'*data;
        
%         options.ReguAlpha = nSmp*options.ReguAlpha;

        switch lower(options.ReguType)
            case {lower('Ridge')}
                for i=1:size(DPrime,1)
                    DPrime(i,i) = DPrime(i,i) + options.ReguAlpha;
                end
            case {lower('Tensor')}
                DPrime = DPrime + options.ReguAlpha*options.regularizerR;
            case {lower('Custom')}
                DPrime = DPrime + options.ReguAlpha*options.regularizerR;
            otherwise
                error('ReguType does not exist!');
        end

        DPrime = max(DPrime,DPrime');
    end
end


[nSmp,nFea] = size(data);

Hb = zeros(nClass,nFea);
for i = 1:nClass,
    index = find(gnd==classLabel(i));
    classMean = mean(data(index,:),1);
    Hb (i,:) = sqrt(length(index))*classMean;
end

elapse.timeW = 0;
elapse.timePCA = cputime - tmp_T;


tmp_T = cputime;

if bPCA
    [dumpVec,eigvalue,eigvector] = svd(Hb,'econ');
    
    eigvalue = diag(eigvalue);
    eigIdx = find(eigvalue < 1e-3);
    eigvalue(eigIdx) = [];
    eigvector(:,eigIdx) = [];

    eigvalue = eigvalue.^2;
    eigvector = eigvector_PCA*(repmat(eigvalue_PCA,1,length(eigvalue)).*eigvector);
else
    WPrime = Hb'*Hb;
    WPrime = max(WPrime,WPrime');

    dimMatrix = size(WPrime,2);
    if Dim > dimMatrix
        Dim = dimMatrix;
    end
    
    if isfield(options,'bEigs')
        if options.bEigs
            bEigs = 1;
        else
            bEigs = 0;
        end
    else
        if (dimMatrix > 1000 & Dim < dimMatrix/10) | (dimMatrix > 500 & Dim < dimMatrix/20) | (dimMatrix > 250 & Dim < dimMatrix/30) 
            bEigs = 1;
        else
            bEigs = 0;
        end
    end
    
    if bEigs
        %disp('use eigs to speed up!');
        option = struct('disp',0);
        if bChol
            option.cholB = 1;
            [eigvector, eigvalue] = eigs(WPrime,R,Dim,'la',option);
        else
            [eigvector, eigvalue] = eigs(WPrime,DPrime,Dim,'la',option);
        end
        eigvalue = diag(eigvalue);
    else
        [eigvector, eigvalue] = eig(WPrime,DPrime);
        eigvalue = diag(eigvalue);

        [junk, index] = sort(-eigvalue);
        eigvalue = eigvalue(index);
        eigvector = eigvector(:,index);

        if Dim < size(eigvector,2)
            eigvector = eigvector(:, 1:Dim);
            eigvalue = eigvalue(1:Dim);
        end
    end
end

for i = 1:size(eigvector,2)
    eigvector(:,i) = eigvector(:,i)./norm(eigvector(:,i));
end


elapse.timeMethod = cputime - tmp_T; 
elapse.timeAll = elapse.timePCA + elapse.timeMethod;
    
