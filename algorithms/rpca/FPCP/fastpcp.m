function [L1, S1, statsPCP] = fastpcp(V, lambda, loops, rank0, rankThreshold, lambdaFactor )

if( nargin < 6 )
  lambdaFactor = 1.0;           
  if( nargin < 5 )
    rankThreshold = 0.01;       
    if( nargin < 4 )
      rank0 = 1;                % initial rank
      if( nargin < 3 )
        loops = 2;              % total number of outer loops
      end
    end
  end
end

if isempty(rankThreshold)
	rankThreshold = 0.01;
end

if isempty(lambdaFactor)
	lambdaFactor = 1.0;
end

% Data size
%[Nrows,Ncols] = size(V);

% Set flag (increments the rank plus one at each iteration)
inc_rank = 1;

% ---------------------------------
% >>>  measure time performance <<<
% ---------------------------------
%t = tic;

% ------------------------
% --- First outer loop ---
  rank = rank0;                     % current rank
  statsPCP.rank(1) = rank;          % save current rank

  % Partial SVD
  %[Ulan Slan Vlan] = lansvd(V, rank, 'L');
  %[Ulan,Slan,Vlan] = svds(V, rank);
  [Ulan,Slan,Vlan] = svdsecon(V, rank); % fastest

  % Current low-rank approximation
  L1 = Ulan*Slan*Vlan';

  % Shrinkage
  S1 = shrink(V-L1, lambda);
  
% ------------------------
% ---    Outer loops   ---

for k = 2:loops,
  
  if(inc_rank == 1)
     lambda = lambda * lambdaFactor;         % modify Lambda at each iteration
     rank = rank + 1;                        % increase rank
  end

  % low rank (partial SVD)
  %[Ulan Slan Vlan] = lansvd(V-S1, rank, 'L');
  %[Ulan,Slan,Vlan] = svds(V-S1, rank);
  [Ulan,Slan,Vlan] = svdsecon(V-S1, rank); % fastest

  currentEvals = diag(Slan);                                    % extract current evals
  statsPCP.rank(k) = length( currentEvals );                    % save current rank
  statsPCP.rho(k) = currentEvals(end) / sum( currentEvals(1:end-1) );       % relative contribution of the last evec

  % simple rule to keep or increase the current rank's value
  if(statsPCP.rho(k) < rankThreshold ) 
     inc_rank = 0;
  else
     inc_rank = 1;
  end

  % Current low-rank approximation
  L1 = Ulan*Slan*Vlan';
  
  % Shrinkage
  S1 = shrink(V-L1, lambda);
end

% ---------------------------------
% >>>  measure time performance <<<
% ---------------------------------
%statsPCP.time = toc(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = shrink(v, lambda)
  u = sign(v).*max(0, abs(v) - lambda);
return 
