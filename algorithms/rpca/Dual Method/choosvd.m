function y = choosvd( n, d)

if n <= 100 
    if d / n <= 0.02
        y = 1;
    else
        y = 0;
    end
elseif n <= 200
    if d / n <= 0.10
        y = 1;
    else
        y = 0;
    end
elseif n <= 300
    if d / n <= 0.13
        y = 1;
    else
        y = 0;
    end
elseif n <= 400
    if d / n <= 0.14
        y = 1;
    else
        y = 0;
    end
elseif n <= 500
    if d / n <= 0.17
        y = 1;
    else
        y = 0;
    end
else
    if d / n <= 0.19
        y = 1;
    else
        y = 0;
    end
end