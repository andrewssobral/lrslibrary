function lambdas = hessianspectrum(problem, x, sqrtprec)
% Returns the eigenvalues of the (preconditioned) Hessian at x.
% 
% function lambdas = hessianspectrum(problem, x)
% function lambdas = hessianspectrum(problem, x, sqrtprecon)
%
% If the Hessian is defined in the problem structure and if no
% preconditioner is defined, returns the eigenvalues of the Hessian
% operator (which needs to be symmetric but not necessarily definite) on
% the tangent space at x. There are problem.M.dim() eigenvalues.
%
% If a preconditioner is defined, the eigenvalues of the composition is
% computed: precon o Hessian. Remember that the preconditioner has to be
% symmetric, positive definite, and is supposed to approximate the inverse
% of the Hessian.
%
% Even though the Hessian and the preconditioner are both symmetric, their
% composition is not symmetric, which can slow down the call to 'eigs'
% substantially. If possible, you may specify the square root of the
% preconditioner as an optional input sqrtprecon. This operator on the
% tangent space at x must also be symmetric, positive definite, and such
% that sqrtprecon o sqrtprecon = precon. Then the spectrum of the symmetric
% operator sqrtprecon o hess o sqrtprecon is computed: it is the same as
% the spectrum of precon o hess, but is generally faster to compute.
% The operator sqrtprecon(x, u[, store]) accepts as input: a point x,
% a tangent vector u and (optional) a store structure.
%
% The input and the output of the Hessian and of the preconditioner are
% projected on the tangent space to avoid undesired contributions of the
% ambient space.
%
% Requires the manifold description in problem.M to have these functions:
% 
%   u_vec = vec(x, u_mat) :
%       Returns a column vector representation of the normal (usually
%       matrix) representation of the tangent vector u_mat. vec must be an
%       isometry between the tangent space (with its Riemannian metric) and
%       a subspace of R^n where n = length(u_vec), with the 2-norm on R^n.
%       In other words: it is an orthogonal projector.
%
%   u_mat = mat(x, u_vec) :
%       The inverse of vec (its adjoint).
%
%   u_mat_clean = tangent(x, u_mat) :
%       Subtracts from the tangent vector u_mat any component that would
%       make it "not really tangent", by projection.
%
%   answer = vecmatareisometries() :
%       Returns true if the linear maps encoded by vec and mat are
%       isometries, false otherwise. It is better if the answer is yes.
%

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, July 3, 2013.
% Contributors: 
% Change log:


    if ~canGetHessian(problem)
        warning('manopt:hessianspectrum:nohessian', ...
                ['The Hessian appears to be unavailable.\n' ...
                 'Will try to use an approximate Hessian instead.\n'...
                 'Since this approximation may not be linear or '...
                 'symmetric,\nthe computation might fail and the '...
                 'results (if any)\nmight make no sense.']);
    end

    vec = @(u_mat) problem.M.vec(x, u_mat);
    mat = @(u_vec) problem.M.mat(x, u_vec);
    tgt = @(u_mat) problem.M.tangent(x, u_mat);
    
    % n: size of a vectorized tangent vector
    % dim: dimension of the tangent space
    % necessarily, n >= dim.
    % The vectorized operators we build below will have at least n - dim
    % zero eigenvalues.
    n = length(vec(problem.M.zerovec(x)));
    dim = problem.M.dim();
    
    % The store structure is not updated by the getHessian call because the
    % eigs function will not take care of it. This might be worked around,
    % but for now we simply obtain the store structure built from calling
    % the cost and gradient at x and pass that one for every Hessian call.
    % This will typically be enough, seen as the Hessian is not supposed to
    % store anything new.
    storedb = struct();
    if canGetGradient(problem)
        [unused1, unused2, storedb] = getCostGrad(problem, x, struct()); %#ok<ASGLU>
    end
    
    hess = @(u_mat) tgt(getHessian(problem, x, tgt(u_mat), storedb));
    hess_vec = @(u_vec) vec(hess(mat(u_vec)));
    
    % Regardless of preconditioning, we can only have a symmetric
    % eigenvalue problem if the vec/mat pair of the manifold is an
    % isometry:
    vec_mat_are_isometries = problem.M.vecmatareisometries();
    
    if ~exist('sqrtprec', 'var') || isempty(sqrtprec)
    
        if ~canGetPrecon(problem)
            
            % There is no preconditinoer : just deal with the (symmetric)
            % Hessian.
            
            eigs_opts.issym = vec_mat_are_isometries;
            eigs_opts.isreal = true;
            lambdas = eigs(hess_vec, n, dim, 'LM', eigs_opts);
            
        else
            
            % There is a preconditioner, but we don't have its square root:
            % deal with the non-symmetric composition prec o hess.
            
            prec = @(u_mat) tgt(getPrecon(problem, x, tgt(u_mat), storedb));
            prec_vec = @(u_vec) vec(prec(mat(u_vec)));
            % prec_inv_vec = @(u_vec) pcg(prec_vec, u_vec);

            eigs_opts.issym = false;
            eigs_opts.isreal = true;
            lambdas = eigs(@(u_vec) prec_vec(hess_vec(u_vec)), ...
                           n, dim, 'LM', eigs_opts);
            
        end
        
    else
        
        % There is a preconditioner, and we have its square root: deal with
        % the symmetric composition sqrtprec o hess o sqrtprec.
        % Need to check also whether sqrtprec uses the store cache or not.

		is_octave = exist('OCTAVE_VERSION', 'builtin');
        if ~is_octave
			narg = nargin(sqrtprec);
		else
			narg = 3;
        end
		
        switch narg
			case 2
				sqrtprec_vec = @(u_vec) vec(tgt(sqrtprec(x, tgt(mat(u_vec)))));
			case 3
				store = getStore(problem, x, storedb);
				sqrtprec_vec = @(u_vec) vec(tgt(sqrtprec(x, tgt(mat(u_vec)), store)));
			otherwise
				error('sqrtprec must accept 2 or 3 inputs: x, u, (optional: store)');
        end
        
        eigs_opts.issym = vec_mat_are_isometries;
        eigs_opts.isreal = true;
        lambdas = eigs(@(u_vec) ...
                      sqrtprec_vec(hess_vec(sqrtprec_vec(u_vec))), ...
                      n, dim, 'LM', eigs_opts);
        
    end

end
