function out = l1proj(Sbar,delta)
sbar=Sbar(:);
n=length(sbar);
norm_sbar=norm(sbar);

if norm_sbar<=delta
    out=0;
else
    flag = 0;
    ssbar = sort(sbar,'descend');
    partialsum=norm_sbar^2;
    for i=1:n-1
        partialsum=partialsum-ssbar(i)^2;
        if i*ssbar(i+1)^2+partialsum<delta^2
            lambda = sqrt((delta^2-partialsum)/i);
            flag = 1;
            break;
        end
    end
    if flag==0
        lambda=delta/sqrt(n);
    end
    out = sign(Sbar).*max(abs(Sbar)-lambda,0);
end
        