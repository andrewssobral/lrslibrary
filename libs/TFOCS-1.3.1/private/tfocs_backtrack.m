% TFOCS_BACKTRACK
% Backtracking helper script.

% Quick exit for no backtracking
if beta >= 1, break; end % SRB changing == to >=

% Quick exit if no progress made
xy = x - y;
xy_sq = tfocs_normsq( xy );
if xy_sq == 0, localL = Inf; break; end

% Compute Lipschitz estimate
if backtrack_simple,
    if isinf( f_x ),
        f_x = apply_smooth( A_x );
    end
    q_x = f_y + tfocs_dot( xy, g_y ) + 0.5 * L * xy_sq;
    localL = L + 2 * max( f_x - q_x, 0 ) / xy_sq;
    backtrack_simple = abs( f_y - f_x ) >= backtrack_tol * max( abs( f_x ), abs( f_y ) );
else
    if isempty( g_Ax ),
        [ f_x, g_Ax ] = apply_smooth( A_x );
    end
    localL = 2 * tfocs_dot( A_x - A_y, g_Ax - g_Ay ) / xy_sq;
end

% Exit if Lipschitz criterion satisfied, or if we hit Lexact
backtrack_steps = backtrack_steps + 1;
if localL <= L || L >= Lexact, break; end
if ~isinf( localL ), 
    L = min( Lexact, localL );
elseif isinf( localL ), localL = L; end
L = min( Lexact, max( localL, L / beta ) );

% TFOCS v1.3 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2013 California Institute of Technology and CVX Research.
% See the file LICENSE for full license information.
