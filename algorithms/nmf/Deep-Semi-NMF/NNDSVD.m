function [W,H] = NNDSVD(A,k,flag);
%
% This function implements the NNDSVD algorithm described in [1] for
% initialization of Nonnegative Matrix Factorization Algorithms.
%
% [W,H] = nndsvd(A,k,flag);
%
% INPUT
% ------------
%
% A    : the input nonnegative m x n matrix A
% k    : the rank of the computed factors W,H
% flag : indicates the variant of the NNDSVD Algorithm
%        flag = 0 --> NNDSVD
%        flag = 1 --> NNDSVDa
%        flag = 2 --> NNDSVDar
%
% OUTPUT
% -------------
%   
% W   : nonnegative m x k matrix
% H   : nonnegative k x n matrix
%
% 
% References:
% 
% [1] C. Boutsidis and E. Gallopoulos, SVD-based initialization: A head
%     start for nonnegative matrix factorization, Pattern Recognition,
%     Elsevier
%
% This code is kindly provided by the authors for research porpuses.
% - Efstratios Gallopoulos (stratis@ceid.upatras.gr)
% - Christos Boutsidis (boutsc@cs.rpi.edu)     
%
% For any problems or questions please send an email to boutsc@cs.rpi.edu 
%--------------------------------------------------------------------------


%----------------------check the input matrix------------------------------
if numel(find(A<0)) > 0
    error('The input matrix contains negative elements !')
end
%--------------------------------------------------------------------------

%size of the input matrix
[m,n] = size(A);

%the matrices of the factorization
W = zeros(m,k);
H = zeros(k,n);

% 1st SVD --> partial SVD rank-k to the input matrix A. 
[U,S,V] = svds(A,k);

%choose the first singular triplet to be nonnegative
W(:,1)     =  sqrt(S(1,1)) * abs(U(:,1) );         
H(1,:)     =  sqrt(S(1,1)) * abs(V(:,1)'); 


% 2nd SVD for the other factors (see table 1 in our paper)
for i=2:k
    uu = U(:,i); vv = V(:,i);
    uup = pos(uu); uun = neg(uu) ;
    vvp = pos(vv); vvn = neg(vv);
    n_uup = norm(uup);
    n_vvp = norm(vvp) ;
    n_uun = norm(uun) ;
    n_vvn = norm(vvn) ;
    termp = n_uup*n_vvp; termn = n_uun*n_vvn;
    if (termp >= termn)
        W(:,i) = sqrt(S(i,i)*termp)*uup/n_uup; 
        H(i,:) = sqrt(S(i,i)*termp)*vvp'/n_vvp;
    else
        W(:,i) = sqrt(S(i,i)*termn)*uun/n_uun; 
        H(i,:) = sqrt(S(i,i)*termn)*vvn'/n_vvn;
    end
end
%------------------------------------------------------------

%actually these numbers are zeros
W((W<0.0000000001))=0.1;
H((H<0.0000000001))=0.1;

% NNDSVDa: fill in the zero elements with the average 
if flag==1
   ind1      =  find(W==0) ;
   ind2      =  find(H==0) ;
   average   =  mean(A(:)) ; 
   W( ind1 ) =  average    ; 
   H( ind2 ) =  average    ;

% NNDSVDar: fill in the zero elements with random values in the space [0:average/100]
elseif flag==2
   ind1      =  find(W==0) ;
   ind2      =  find(H==0) ;
   n1        =  numel(ind1);
   n2        =  numel(ind2);
   
   average   =  mean(A(:))       ;
   W( ind1 ) =  (average*rand(n1,1)./100)  ; 
   H( ind2 ) =  (average*rand(n2,1)./100)  ;   
end

end
%--------------------------------------------------------------------------
%end of the nndsvd function



%This function sets to zero the negative elements of a matrix
%--------------------------------------------------------------------------
function [Ap] = pos(A)
Ap = (A>=0).*A;
end
%--------------------------------------------------------------------------

%This functions sets to zero the positive elements of a matrix and takes
%the absolute value of the negative elements
%--------------------------------------------------------------------------
function [Am] = neg(A);
Am = (A<0).*(-A);
end
%--------------------------------------------------------------------------
