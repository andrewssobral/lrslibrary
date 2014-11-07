%% [matrix matrix] tfocs_interface(matrix)
% M - 2d matrix
% L - low rank matrix
% S - sparse matrix
%
function [L, S] = tfocs_interface(M,constraint)
  X = M;
  %lambda = 1/sqrt(min(size(M,1),size(M,2))); 
  lambda  = 1e-2;
  largescale = false;
  opts = [];
  opts.stopCrit       = 4;
  opts.printEvery     = 1;
  opts.tol            = 1e-4;
  opts.maxIts         = 25;
  opts.errFcn{1}      = @(f,d,p) norm(p{1}+p{2}-X,'fro')/norm(X,'fro');
  
  for inequality_constraints = 0:1
      if inequality_constraints
          % if we already have equality constraint solution,
          %   it would make sense to "warm-start":
  %         x0      = { LL_0, SS_0 };
          % but it's more fair to start all over:
          x0      = { X, zeros(size(X))   };
          z0      = [];
      else
          x0      = { X, zeros(size(X))   };
          z0      = [];
      end


      obj    = { prox_nuclear(1,largescale), prox_l1(lambda) };
      affine = { 1, 1, -X };

      mu = 1e-4;
      if inequality_constraints
          epsilon  = 0.5;
          dualProx = prox_l1(epsilon);
      else
          dualProx = proj_Rn;
      end

      %tic
      % call the TFOCS solver:
      [x,out,optsOut] = tfocs_SCD(obj, affine, dualProx, mu, x0, z0, opts);
      %toc

      % save the variables
      LL = x{1};
      SS = x{2};
      if ~inequality_constraints
          z0      = out.dual;
          LL_0    = LL;
          SS_0    = SS;
      end

  end % end loop over "inequality_constriants" variable

  L = [];
  S = [];
  
  % equality constraints
  if(constraint == 1)
    L = LL_0;
    S = SS_0;
  end
  
  % inequality constraints
  if(constraint == 2)
    L = LL;
    S = SS;
  end
end

