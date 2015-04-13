function [tracedtensor] = diagsum(tensor1, d1, d2)
% C = DIAGSUM(A, d1, d2) Performs the trace
% C(i[1],...,i[d1-1],i[d1+1],...,i[d2-1],i[d2+1],...i[n]) =
%              A(i[1],...,i[d1-1],k,i[d1+1],...,i[d2-1],k,i[d2+1],...,i[n])
% (Sum on k).
%
% C = DIAGSUM(A, d1, d2) traces A along the diagonal formed by dimensions d1
% and d2. If the lengths of these dimensions are not equal, DIAGSUM traces
% until the end of the shortest of dimensions d1 and d2 is reached. This is
% an analogue of the built in TRACE function.
%
% Wynton Moore, January 2006


dim1=size(tensor1);
numdims=length(dim1);


%check inputs
if d1==d2
    tracedtensor=squeeze(sum(tensor1,d1));
elseif numdims==2
    tracedtensor=trace(tensor1);
elseif dim1(d1)==1 && dim1(d2)==1
    tracedtensor=squeeze(tensor1);
else


    %determine correct permutation
    swapd1=d1;swapd2=d2;
    
    if d1~=numdims-1 && d1~=numdims && d2~=numdims-1
        swapd1=numdims-1;
    elseif d1~=numdims-1 && d1~=numdims && d2~=numdims
        swapd1=numdims;
    end
    if d2~=numdims-1 && d2~=numdims && swapd1~=numdims-1
        swapd2=numdims-1;
    elseif d2~=numdims-1 && d2~=numdims && swapd1~=numdims
        swapd2=numdims;
    end
    
    
    %prepare for construction of selector tensor
    temp1=eye(numdims);
    permmatrix=temp1;
    permmatrix(:,d1)=temp1(:,swapd1);
    permmatrix(:,swapd1)=temp1(:,d1);
    permmatrix(:,d2)=temp1(:,swapd2);
    permmatrix(:,swapd2)=temp1(:,d2);

    selectordim=dim1*permmatrix;
    permvector=(1:numdims)*permmatrix;


    %construct selector tensor
    if numdims>3
        selector = ipermute(outer(ones(selectordim(1:numdims-2)), ...
                                  eye(selectordim(numdims-1), ...
                                      selectordim(numdims)), ...
                                  0), ...
                            permvector);
    else
        %when numdims=3, the above line gives ndims(selector)=4. This
        %routine avoids that error. When used with GMDMP, numdims will be
        %at least 4, so this routine will be unnecessary.
        selector2=eye(selectordim(numdims-1), selectordim(numdims));
        selector=zeros(selectordim);
        for j=1:selectordim(1)
            selector(j, :, :)=selector2;
        end
        selector=ipermute(selector, permvector);
    end
    
    
    %perform trace, discard resulting singleton dimensions
    tracedtensor=sum(sum(tensor1.*selector, d1), d2);
    tracedtensor=squeeze(tracedtensor);
	
    
end


%correction for abberation in squeeze function:
%size(squeeze(rand(1,1,2)))=[2 1]
nontracedimensions=dim1;
nontracedimensions(d1)=[];
if d2>d1
    nontracedimensions(d2-1)=[];
else
    nontracedimensions(d2)=[];
end
tracedsize=size(tracedtensor);
% Next line modified, Nicolas Boumal, April 30, 2012, such that diagsum(A,
% 1, 2) would compute the trace of A, a 2D matrix.
if length(tracedsize)==2 && tracedsize(2)==1 && ...
   (isempty(nontracedimensions) || tracedsize(1)~=nontracedimensions(1))

    tracedtensor=tracedtensor.';
    
end
