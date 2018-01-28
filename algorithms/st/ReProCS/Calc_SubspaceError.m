function[SE_err] = Calc_SubspaceError(A, B)
%% Calculate subspace error of two subspaces

%Calculate the basis (orthonormal column vectors) of span(A), span(B)
if(~isempty(A) && ~isempty(B))
    [QA, ~] = qr(A);
    [QB, ~] = qr(B);
    
    ra = rank(A);
    rb = rank(B);
    
    [m,~] = size(A);
    SE_err = norm((eye(m) - QA(:, 1 : ra)  * QA(:, 1 : ra)') ...
        * QB(:, 1 : rb), 2);
else
    SE_err = eps;
end
end