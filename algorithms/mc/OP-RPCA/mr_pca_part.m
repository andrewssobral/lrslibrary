function [L,C] = mr_pca_part(M,Omega,lambda, varargin)
% This is a function to solve mr_PCA with partial observation
% M is the input matrix with partical observation
% Omega is the sampling set, i.e., Omega_(i,j)=1 if M_i,j is observed
% lambda is the input paremter
% L is the output low rank matrix
% C is the output column sparse matrix
% 'full_svd' if 1, use a full svd, otherwise a partial one
%
[m, n]=size(M);
if size(Omega)~=size(M)
    error('Size of Omega and M are incosistent');
end    
rank_L=ceil(min(m, n)*0.1);
increaseK=10;
f_svd=1;
if (mod(length(varargin), 2) ~= 0 ),
    error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
end
parameterCount = length(varargin)/2;
for parameterIndex = 1:parameterCount,
    parameterName = varargin{parameterIndex*2 - 1};
    parameterValue = varargin{parameterIndex*2};
    switch lower(parameterName)
        case 'full_svd'
            f_svd=parameterValue;
        case 'increaseK'
            increaseK=parameterValue;
        case 'ini_rank'
            rank_L=parameterValue;
 otherwise
            error(['Sorry, the parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
    end
end

delta=0.00001;
mu_temp=0.99*norm(M);
mu_bar=delta*mu_temp;
eta=0.9;
tol=0.000001*norm(M,'fro');
stopping=0;

% L_temp0, C_temp0 is L_k and C_k
% L_temp1, C_temp1 is L_{k-1} and C_{k-1}
L_temp0=zeros(m,n);
L_temp1=zeros(m,n); 
C_temp0=zeros(m,n);
C_temp1=zeros(m,n);
t_temp0=1;
t_temp1=1;
k=0;

while stopping~=1
  YL=L_temp0+(t_temp1-1)/t_temp0*(L_temp0-L_temp1);
  YC=C_temp0+(t_temp1-1)/t_temp0*(C_temp0-C_temp1);
  
  M_difference=(YL+YC-M).*Omega;

  GL=YL-0.5*(M_difference);
  [L_new, rank_L]=iterate_L(GL, mu_temp/2, f_svd, increaseK, rank_L+1);

  GC=YC-0.5*(M_difference);
  C_new=iterate_C(GC, mu_temp*lambda/2);

  t_new=(1+sqrt(4*t_temp0^2+1))/2;
  mu_new=max(eta*mu_temp, mu_bar);
  
  %%%%%%% Now to decide whether to stop
  S_L=2*(YL-L_new)+(L_new+C_new-YL-YC);
  S_C=2*(YC-C_new)+(L_new+C_new-YL-YC);
  if norm(S_L, 'fro')^2+norm(S_C, 'fro')^2 <= tol^2;
      stopping=1;
  else
    L_temp1=L_temp0;
    L_temp0=L_new;
    C_temp1=C_temp0;
    C_temp0=C_new;
    t_temp1=t_temp0;
    t_temp0=t_new;
    mu_temp=mu_new;
    k=k+1
  end  
end
L=L_new;
C=C_new;

%%%%%%%%%%%%%%%%%%%%%%%%%
function [output,rank_out]=iterate_L(L, epsilon, f_svd, increaseK, starting_K) 
if f_svd~=1
    a.minSingValue=epsilon;
    a.increaseK=increaseK;
    [U, S, V]=lansvd(L, starting_K, 'T', a);
    rank_out=min(size(S));
else
   [U, S, V]=svd(L);
   rank_out=0;
end
   for i=1:min(size(S))
   if S(i,i)> epsilon
       S(i,i)=S(i,i)-epsilon;
   elseif S(i,i)<-epsilon
       S(i,i)=S(i,i)+epsilon;
   else 
       S(i,i)=0;
   end
end   
output=U*S*V';
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output=iterate_C(C, epsilon)
[m,n]=size(C);
for i=1:n
    temp=C(:,i);
    norm_temp=norm(temp);
    if norm_temp>epsilon
        temp=temp-temp*epsilon/norm_temp;
    else
        temp=zeros(m,1);
    end
    output(:,i)=temp;
end
return;

