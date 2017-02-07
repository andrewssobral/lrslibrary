%%
% [Idx, Omega] = subsampling(data, obs)
%
% data - matrix/tensor
% obs - double (percentual of observed entries) [0...1]
%
% Idx - observed indexes
% Omega - binary matrix/tensor
%
function [Idx, Omega] = subsampling(data, obs)
  if(nargin < 1 || nargin > 2)
    error('data must be defined');
  end
  if(nargin == 1)
    obs = 0.5; % 50% of observed values (default)
  end
  nae = numel(double(data));
  rp = randperm(nae)';
  k = floor(obs*nae);
  Idx = rp(1:k);
  Omega = zeros(size(data));
  Omega(Idx) = 1;
end
