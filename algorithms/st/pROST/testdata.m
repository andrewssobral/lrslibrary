% Copyright (c) Technische Universität München, 2012
% All rights reserved.
%
% Author: Clemens Hage
% Contact: hage@tum.de
% Date: 2012-11-12
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

function [X,L,S]=testdata(m,n,k,rho,amplitude,sigma)
% This function creates a random data set containing of n samples of
% m-dimensional data. The data is centered and normally distributed with a
% sample standard deviation of 1.

% SVD is applied to the dataset and all singular values sigma_i with i>k
% are set to zero, thus establishing a k-dimensional subspace.

% In order to perturb the dataset with sparse outliers, the variable rho can
% be set to a value between 0 and 1 which corresponds to the density of
% the uniformly distributed entries. The outliers are uniformly distributed
% in [-amplitude amplitude]

% in case you want to add Gaussian noise, choose non-zero value for sigma.

% start off with random mxn matrix
R=randn(m,n);

[U,Sigma,V]=svd(R,0);

% set singular values k+1 to m to zero
for i=k+1:m
     Sigma(i,i)=0;
end

% low-rank data matrix
L=U*Sigma*V;

% normalize to sample variance of 1
L=L./std(L(:));

% Gaussian noise
N = sigma*randn(size(L));

% sparse outliers
S=mysprand(m,n,rho,amplitude);

X = L + N + S;

end

function S = mysprand(m,n,rho,amplitude)
% this does essentially the same as MATLAB's default function sprandn, just
% that the number of nonzero entries is precise in contrast to MATLAB's.
% Also, you can control the amplitude of the entries. Here we choose a
% zero-mean uniform distribution between [-amplitude amplitude]


% total number of entries
elements = m * n;

% number of nonzero entries
nzentries = ceil(rho * m * n);

% generate the nonzero entries
values = amplitude * 2 * (rand(nzentries,1)-0.5);

% random positions
indices = randperm(elements,nzentries);

% convert to indices
[i,j]=ind2sub([m n],indices);

% create sparse matrix
S = sparse(i,j,values,m,n,nzentries);

end