% Author Bruno Luong <brunoluong@yahoo.com>
% Last update: 07/April/2009
% Script to tets min/max selection

clear

try
    minkmex(1,1);
    minkmex(1,1);
catch
    minmax_install();
end

n=1e7;
k=10;

ntest=10;

% Timing
disp('Time the algorithms for few seconds...');
tmink=zeros(1,ntest);
tmaxk=zeros(1,ntest);
tsort=zeros(1,ntest);
for i=1:ntest
    list=rand(1,n);
    
    tic
    mn=mink(list,k);
    tmink(i)=toc;
    
    tic
    mx=maxk(list,k);
    tmaxk(i)=toc;
    
    tic
    s=sort(list);
    smn=s(1:k);
    smx=s(end:-1:end-k+1);
    tsort(i)=toc;
    
    if ~isequal(mn,smn) || ~isequal(mx,smx)
        keyboard;
    end
end

fprintf('Timing mink: %f [s]\n',median(tmink));
fprintf('Timing maxk: %f [s]\n',median(tmaxk));
fprintf('Timing sort: %f [s]\n',median(tsort));