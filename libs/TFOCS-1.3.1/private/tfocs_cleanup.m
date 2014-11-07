%SOLVER_CLEANUP    TFOCS helper script
%   Performs the final set of operations for our templated solvers: performs
%   final calculations, prints final outputs, cleans up unused data.

if saddle,
    if isempty(g_Ax),
        [ f_x, g_Ax ] = apply_smooth(A_x);
    end
    out.dual = get_dual( g_Ax );
    % Oct 13: make it compatible:
    if isa( out.dual, 'tfocs_tuple')
        out.dual = cell( out.dual );
    end
elseif isinf(f_x),
    f_x = apply_smooth(A_x);
end
if isinf( C_x ),
    C_x = apply_projector( x );
end    
if fid && printEvery,
	fprintf( fid, 'Finished: %s\n', status );
end
out.niter = n_iter;
out.status = status;
d.niter = 'Number of iterations';
if saveHist,
    out.f(n_iter) = maxmin * ( f_x + C_x );
    out.f(n_iter+1:end) = [];
    out.normGrad(n_iter+1:end) = [];
    out.stepsize(n_iter+1:end) = [];
    out.theta(n_iter+1:end,:) = [];
    if countOps,
        out.counts(n_iter+1:end,:) = [];
    end
    if ~isempty(errFcn),
        out.err(n_iter+1:end,:) = [];
    end
    d.f        = 'Objective function history';
    d.normDecr = 'Decrement norm';
    d.stepsize = 'Stepsize';
    d.theta    = 'Acceleration parameter history';
    if countOps,
        d.counts   = strvcat(...
        'k x 4 arry, with columns [F,G,A,P] where',...
        'F: Number of function evaluations of the smooth function',...
        'G: Number of gradient evaluations of the smooth function',...
        'A: Number of calls to the linear operator and its transpose',...
        'N: Number of calls to the nonsmooth function (w/o projection)',...
        'P: Number of calls to the projection operator' );
    end
    if ~isempty(errFcn)
        d.err = 'Error, determined by evaluating the user-supplied error function';
    end
end
out.description = d;
if countOps,
    clear global tfocs_count___
end
    
% TFOCS v1.3 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2013 California Institute of Technology and CVX Research.
% See the file LICENSE for full license information.
