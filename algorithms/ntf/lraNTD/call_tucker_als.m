function [Y out]=call_tucker_als(X,opts)
% call TUCKER_ALS in tensor toolbox 2.4 in TDALAB
% Outputs: 
%           Y: ktensor
%           A: loading factors
%           G: tucker core
%
%TUCKER_ALS Higher-order orthogonal iteration.
%
%   T = TUCKER_ALS(X,R) computes the best rank(R1,R2,..,Rn)
%   approximation of tensor X, according to the specified dimensions
%   in vector R.  The input X can be a tensor, sptensor, ktensor, or
%   ttensor.  The result returned in T is a ttensor.
%
%   T = TUCKER_ALS(X,R,'param',value,...) specifies optional parameters and
%   values. Valid parameters and their default values are:
%      'tol' - Tolerance on difference in fit {1.0e-4}
%      'maxiters' - Maximum number of iterations {50}
%      'dimorder' - Order to loop through dimensions {1:ndims(A)}
%      'init' - Initial guess [{'random'}|'nvecs'|cell array]
%      'printitn' - Print fit every n iterations {1}
%
%   [T,U0] = TUCKER_ALS(...) also returns the initial guess.
%
%   Examples:
%   X = sptenrand([5 4 3], 10);
%   T = tucker_als(X,2);        %<-- best rank(2,2,2) approximation 
%   T = tucker_als(X,[2 2 1]);  %<-- best rank(2,2,1) approximation 
%   T = tucker_als(X,2,'dimorder',[3 2 1]);
%   T = tucker_als(X,2,'dimorder',[3 2 1],'init','nvecs');
%   U0 = {rand(5,2),rand(4,2),[]}; %<-- Initial guess for factors of T
%   T = tucker_als(X,2,'dimorder',[3 2 1],'init',U0);
%
%   See also TTENSOR, TENSOR, SPTENSOR, KTENSOR.
%
%MATLAB Tensor Toolbox.
%Copyright 2010, Sandia Corporation. 

% This is the MATLAB Tensor Toolbox by Brett Bader and Tamara Kolda. 
% http://csmr.ca.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2010) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in tensor_toolbox/LICENSE.txt
% $Id: tucker_als.m,v 1.8 2010/03/19 23:46:32 tgkolda Exp $

N=ndims(X);
if ~exist('opts','var')
    opts = struct;
end
defoptions = struct('NumOfComp',[],'tol',1e-12,'maxiters',50,...
    'init','random','printitn',10,'dimorder',1:N);
if ~exist('opts','var')
    opts = struct;
end
[R,tol,maxiters,init,printitn,dimorder] = scanparam(defoptions,opts);
if strcmpi(init,'load')
    init=load('TDinit.mat','Ainit');
    init=Ainit;
    clear Ainit;
    if isempty(init)
        fprintf('[TDALAB] CP_init.mat is empty. Please check the ''init'' parameter.\n');
        return;
    end
end
% save('Ynoi.mat','X','-v7.3');
[Y]=tucker_als(X,R,'init',init,'maxiters',maxiters,'tol',tol,'dimorder',dimorder,...
    'printitn',printitn);