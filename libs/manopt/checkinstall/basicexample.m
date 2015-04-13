function basicexample
    
    % Verify that Manopt was indeed added to the Matlab path.
    if isempty(which('spherefactory'))
        error(['You should first add Manopt to the Matlab path.\n' ...
		       'Please run importmanopt first.']);
    end
    
    % Generate the problem data.
    n = 1000;
    A = randn(n);
    A = .5*(A+A');
    
    % Create the problem structure.
    manifold = spherefactory(n);
    problem.M = manifold;
    
    % Define the problem cost function and its gradient.
    problem.cost = @(x) -x'*(A*x);
    problem.grad = @(x) manifold.egrad2rgrad(x, -2*A*x);
    
    % Numerically check gradient consistency.
    checkgradient(problem);
 
    % Solve.
    % The trust-regions algorithm requires the Hessian. Since we do not
    % provide it, it will go for a standard approximation of it. The first
    % instruction tells Manopt not to issue a warning when this happens.
    warning('off', 'manopt:getHessian:approx');
    [x xcost info] = trustregions(problem); %#ok<ASGLU>
    
    % Display some statistics.
    figure;
    semilogy([info.iter], [info.gradnorm], '.-');
    xlabel('Iteration #');
    ylabel('Gradient norm');
    title('Convergence of the trust-regions algorithm on the sphere');
    
end
