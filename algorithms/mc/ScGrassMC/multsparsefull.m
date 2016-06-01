function A = multsparsefull(C, B, mask)
% FUNCTION A = MULTSPARSEFULL(C, B, MASK)
%
% Efficient multiplication of a sparse and a full matrix. This code relies
% on an ugly C-mex trick, and should not be used for anything else than
% what it was written for.
%
% D is a sparse m-by-n matrix with entries C at positions specified by the
%        sparse mask matrix. Please note that the entries in mask *will*
%        be modified by this function: it is a dummy place holder.
% B is a full n-by-p matrix.
% A is the full m-by-p matrix DB.
% Complexity: O( p*nnz(D) ).
%
% Nicolas Boumal, UCLouvain, May 20, 2011.
% http://perso.uclouvain.be/nicolas.boumal/RTRMC/
%
% SEE ALSO: multfullsparse

    setsparseentries(mask, C)
    A = mask*B;

end
