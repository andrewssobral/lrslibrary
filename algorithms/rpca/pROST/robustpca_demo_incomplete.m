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

% This script runs a quick numerical demo on Robust PCA (batch mode with
% incomplete observation)

% Dimension of the data set
m = 400;

% Number of samples
n = 400;

% Dimension of the underlying subspace
k = 40; % i.e. relative rank of 0.1

% Relative density of outliers
rho = 0.1; % i.e. 10% of the matrix entries are corrupted

% Max. magnitude of the outliers
amplitude = 5;

% Create test data
[X,L_0]=testdata(m,n,k,rho,amplitude,0);

support=0.5; % 50% of the entries are observed

% Create support set

% number of observed entries
nentries = ceil(support * m * n);

% random positions
indices = randperm(m*n,nentries);

% convert to indices
[i,j]=ind2sub([m n],indices);

% create sparse matrix
Omega = sparse(i,j,true(nentries,1),m,n,nentries);

% run Robust PCA
L=robustpca_batch(X,k,'atansquare','L_0',L_0,'Omega',Omega);

% Compute relative subspace reconstruction error
err = norm(Omega.*L-Omega.*L_0,'fro') / norm(Omega.*L_0,'fro');
display(['Relative subspace reconstruction error (robust method): ' num2str(err)])

% % Compare to common PCA approach
X_sampled = Omega.*X;
[U,~,~]=svds(X_sampled,k);

L_PCA = U * (U' * X_sampled);

err_PCA = norm(Omega.*L_PCA-Omega.*L_0,'fro') / norm(Omega.*L_0,'fro');
display(['Relative subspace reconstruction error (common PCA): ' num2str(err_PCA)])



